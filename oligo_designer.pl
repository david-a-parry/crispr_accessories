#!/usr/bin/perl
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

#TO DO - CHECK WHETHER PAM SITE IS INTRONIC
#ALLOW FOR INTRONIC MUTATIONS
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use DBI;
use LWP::Simple;
use XML::Simple;
use Bio::SeqIO;
use Bio::DB::GenBank;

my $gb;     #input genbank file
my $pos;    #position in seq to mutate
my $int_pos = 0;    #intron/UTR position if not coding
my $mut;            #sequence to change
my $flanks    = 60; #flanks for repair template
my $min_score = 50;
my $transcript;     #target transcript RefSeq id
my $f_prefix = "CACC"; #5' restriction site compatible addition to forward guide
my $r_prefix = "AAAC"; #5' restriction site compatible addition to rev guide
my $extra_bases = "G"
  ; #extra bases to add after prefix of forward guide and at 3' end of rev guide
my $splice_distance = 8; #print a warning if PAM site mutation in repair template is within this many bp of a splice site
my $build  = 'hg19';
my %codons = (
    AAA => "K",
    AAC => "N",
    AAG => "K",
    AAT => "N",
    ACA => "T",
    ACC => "T",
    ACG => "T",
    ACT => "T",
    AGA => "R",
    AGC => "S",
    AGG => "R",
    AGT => "S",
    ATA => "I",
    ATC => "I",
    ATG => "M",
    ATT => "I",
    CAA => "Q",
    CAC => "H",
    CAG => "Q",
    CAT => "H",
    CCA => "P",
    CCC => "P",
    CCG => "P",
    CCT => "P",
    CGA => "R",
    CGC => "R",
    CGG => "R",
    CGT => "R",
    CTA => "L",
    CTC => "L",
    CTG => "L",
    CTT => "L",
    GAA => "E",
    GAC => "D",
    GAG => "E",
    GAT => "D",
    GCA => "A",
    GCC => "A",
    GCG => "A",
    GCT => "A",
    GGA => "G",
    GGC => "G",
    GGG => "G",
    GGT => "G",
    GTA => "V",
    GTC => "V",
    GTG => "V",
    GTT => "V",
    TAA => "*",
    TAC => "Y",
    TAG => "*",
    TAT => "Y",
    TCA => "S",
    TCC => "S",
    TCG => "S",
    TCT => "S",
    TGA => "*",
    TGC => "C",
    TGG => "W",
    TGT => "C",
    TTA => "L",
    TTC => "F",
    TTG => "L",
    TTT => "F"
);
my $help;

GetOptions(
    "genbank=s"          => \$gb,
    "transcript=s"       => \$transcript,
    "coordinate=i"       => \$pos,
    "intron_position=i"  => \$int_pos,
    "mutation=s"         => \$mut,
    "f|flanks=i"         => \$flanks,
    "score=i"            => \$min_score,
    "p|forward_prefix=s" => \$f_prefix,
    "reverse_prefix=s"   => \$r_prefix,
    "extra=s"            => \$extra_bases,
    "distance=i"         => \$splice_distance,
    "build=s"            => \$build,
    "help"               => \$help,
) or usage("Error retriveing user options.");

usage()                                       if $help;
usage("--genbank argument is required.\n")    if not $gb;
usage("--transcript argument is required.\n") if not $transcript;
usage("--coordinate argument is required.\n") if not $pos;
usage("--mutation argument is required.\n")   if not $mut;

check_mutation($mut);
my %transcript = get_gene_coordinates($transcript);
my $pos_coord = get_cdna_genomic_coordinate( \%transcript, $pos, $int_pos );
print "$transcript{id}: c.$pos"; 
if ($int_pos){
    if ($int_pos > 0){
        print "+$int_pos";
    }else{
        print "$int_pos";
    }
}
print " matches genomic coordinate $transcript{chrom}:$pos_coord\n";
my $stream_obj = Bio::SeqIO->new(
    -file   => "<$gb",
    -format => "genbank",
);

my %guides = ();
my $closest_guide;
my $closest;
my $pam_coord;
my $target_is_on_f_strand = 0;
while ( my $seq_object = $stream_obj->next_seq ) {
    my $dna = $seq_object->seq;
    my @matches = search_sequence( $transcript{dna}, $dna );
    $target_is_on_f_strand++;
    push @matches, search_sequence( $transcript{dna}, revcomp($dna) );
    if ( @matches < 1 ) {
        die
"ERROR - No matches found for target sequence in $transcript{id} - please check your sequence!\n";
    }
    if ( @matches > 1 ) {
        my $m = scalar @matches;
        die
"ERROR - Multiple ($m) matches found for target sequence in $transcript{id} - please use a specific sequence!\n";
    }
    my $target_coord =
      $transcript{starts}->[0] + ( $matches[0] - ( length($dna) ) ) - 1;
    print "Genbank target sequence matches $transcript{chrom}:$target_coord-"
      . ( $target_coord + length($dna) - 1 ) . "\n";

    foreach my $guide ( grep { $_->primary_tag eq 'protein_bind' }
        $seq_object->get_SeqFeatures )
    {
        my ( $strand, $id, $seq, $score );

        foreach my $tag ( $guide->get_all_tags ) {
            if ( $tag eq 'bound_moiety' ) {
                ( $strand, $id ) = parse_moiety( $guide->get_tag_values($tag) );
                if ( not defined $strand or not defined $id ) {
                    die
"Error parsing guide tag $tag, value: $guide->get_tag_values($tag)\n";
                }
            }
            elsif ( $tag eq 'note' ) {
                ( $score, $seq ) = parse_note( $guide->get_tag_values($tag) );
                if ( not defined $score or not defined $seq ) {
                    die
"Error parsing guide tag $tag, value: $guide->get_tag_values($tag)\n";
                }
            }
        }
        next if $score < $min_score;
        $guides{$id}->{seq}    = $seq;
        $guides{$id}->{strand} = $strand;
        $guides{$id}->{score}  = $score;
    }

    if ( not keys %guides ) {
        die "No guides have a score of $min_score or higher\n";
    }

    #find closest remaining guide RNA
    foreach my $id ( sort { $guides{$b}->{score} <=> $guides{$a}->{score} }
        keys %guides )
    {
        my @found = ();
        my $guide_is_on_forward_strand;
        @found = ( search_sequence( $transcript{dna}, $guides{$id}->{seq} ) );
        $guide_is_on_forward_strand++ if @found;
        push @found,
          ( search_sequence( $transcript{dna}, revcomp( $guides{$id}->{seq} ) )
          );
        if ( @found < 1 ) {
            die
"ERROR - Guide RNA $id does not match target sequence - please check your sequence!\n";
        }
        if ( @found > 1 ) {
            my $m = scalar @found;
            die
"ERROR - Multiple ($m) matches found for Guide RNA $id - please check your sequences!\n";
        }
        my $g_start =
          $found[0] -
          ( length( $guides{$id}->{seq} ) ) +
          $transcript{starts}->[0] - 1;
        my $g_end = $found[0] + $transcript{starts}->[0] - 2;
        print "Guide $id matches $transcript{chrom}:$g_start-$g_end\n";
        if ( !$guide_is_on_forward_strand ) {
            $guides{$id}->{pam_coord}                  = $g_start;
            $guides{$id}->{guide_is_on_reverse_strand} = 1;
        }
        else {
            $guides{$id}->{pam_coord}                  = $g_end;
            $guides{$id}->{guide_is_on_reverse_strand} = 0;
        }
        if ( not defined $closest_guide ) {
            $closest_guide = $id;
            $closest       = abs( $guides{$id}->{pam_coord} - $pos_coord );
            $guides{$id}->{pam_dist} = $closest;
        }
        else {
            my $dist = abs( $guides{$id}->{pam_coord} - $pos_coord );
            $guides{$id}->{pam_dist} = $dist;
            if ( $closest > $dist ) {
                $closest_guide = $id;
                $closest       = $dist;
            }
        }
    }
    print
"\nClosest guide is $closest_guide (PAM site is $closest from coordinate ($pos_coord) at $guides{$closest_guide}->{pam_coord})\n\n";
}

foreach my $id ( sort { $guides{$b}->{pam_dist} <=> $guides{$a}->{pam_dist} }
    keys %guides )
{
    print "Guide $id:\nPAM distance is $guides{$id}->{pam_dist} from "
      . " $pos_coord at $guides{$id}->{pam_coord}\n";
    print "Score = $guides{$id}->{score}\n";
    print "Forward Guide: $f_prefix$extra_bases$guides{$id}->{seq}\n";
    print "Reverse Guide: $r_prefix"
      . ( revcomp( $guides{$id}->{seq} ) )
      . ( revcomp($extra_bases) ) . "\n";

    #print Dumper %guides;
    my $template =
      get_repair_template( \%transcript, $guides{$id}, $pos_coord, $mut );
    if ($template) {
        print "Repair template: $template\n";
    }
    else {
        print "No synonymous PAM site mutation found for $id.\n\n";
    }
    print "\n";
}

############################################
sub check_mutation {
    my $m = shift;
    if ( $m =~ /[^acgtACGT]/ ) {
        die "Non-DNA characters detected in mutation $m\n";
    }
    if ( length($m) != 1 ) {
        die
"Only SNVs are currently supported. Mutation must be a single letter.\n";
    }
}
############################################
sub get_repair_template {
    my ( $t, $guide, $pos, $mut ) = @_;
    my $mod_dna   = $t->{dna};
    my $dna_start = $t->{starts}->[0];

    #get region including flanks spanning PAM site and mutation position
    my $min = $guide->{pam_coord} < $pos ? $guide->{pam_coord} : $pos;
    my $max = $guide->{pam_coord} > $pos ? $guide->{pam_coord} : $pos;
    $min -= $flanks;
    $max += $flanks;

    #GET COORDINATES OF THE 1st G of the PAM SITE
    my $first_g_coord;

    #the template sequence in genbank file could be in either orientation
    #need to determine our guides orientation relative to genomic sequence
    if ( $guide->{guide_is_on_reverse_strand} ) {
        $first_g_coord = $guide->{pam_coord};
    }
    else {
        $first_g_coord = $guide->{pam_coord} - 1;
    }

    #change $mod_dna to our mutant version. Won't work for deletions
    my $changed_base = substr( $mod_dna, $pos - $dna_start, 1, uc($mut) );
    if ( $changed_base eq $mut ) {
        die "ERROR The mutated residue ($mut) is the same as reference!\n";
    }

    #CHECK our mutation hasn't already removed guide site
    if ( $first_g_coord == $pos or $first_g_coord + 1 == $pos ) {
        print
"Your mutation removes PAM site for this guide; no further mutation necessary.\n";
        return substr( $mod_dna, $min - $dna_start, $max - $min );
    }

    #CHECK CODONS FOR EACH G
    foreach my $g ( $first_g_coord, $first_g_coord + 1 ) {

        #need to pass $mod_dna in case our mutation has changed
        #pam site

        #CHECK WHETHER SITE IS INTRONIC/UTR
        if ( !is_coding( $t, $g ) ) {
            if ( is_exonic( $t, $g ) ) {
                my $utr_pos;
                if ( $t->{cdsStart} > $g ) {
                    if ( $t->{strand} eq '+' ) {
                        $utr_pos = "c." . ( $g - $t->{cdsStart} );
                    }
                    else {
                        $utr_pos = "c.*" . ( abs( $g - $t->{cdsStart} ) );
                    }
                }
                elsif ( $t->{cdsEnd} < $g ) {
                    if ( $t->{strand} eq '+' ) {
                        $utr_pos = "c.*" . ( $g - $t->{cdsEnd} );
                    }
                    else {
                        $utr_pos = "c.-" . ( $g - $t->{cdsEnd} );
                    }
                }
                else {
                    die "Couldn't determine UTR position for $t->{chrom}:$g\n";
                }
                my $m = 'A';    #arbitrarily change to an A
                my $nt = substr( $mod_dna, $g - $dna_start, 1, uc($m) );
                die "ERROR - PAM position is not a G/C! "
                  if lc($nt) ne 'g' and lc($nt) ne 'c';
                print "Changing $t->{chrom}:$g from $nt to $m (c.$utr_pos)\n";
                return substr( $mod_dna, $min - $dna_start, $max - $min );
            }
            elsif ( is_intronic( $t, $g ) ) {
                my ( $c_pos, $i_pos ) = get_intronic_position( $t, $g );
                my $coord_string = "c.$c_pos";
                if ($i_pos > 0){
                    $coord_string .= "+$i_pos";
                }else{
                    $coord_string .= "$i_pos";
                }
                my $m = 'A';    #arbitrarily change to an A
                my $nt = uc(substr( $mod_dna, $g - $dna_start, 1, uc($m) ));
                die "ERROR - PAM position is not a G/C! "
                  if lc($nt) ne 'g' and lc($nt) ne 'c';
                if (abs($i_pos) <= $splice_distance){
                    print "WARNING - intronic PAM site mutation is only " . abs($i_pos) . " away from a splice site!\n";
                }
                print
"Changing $t->{chrom}:$g from $nt to $m ($coord_string$nt>$m) for intronic PAM site mutation of repair template.\n";
                return substr( $mod_dna, $min - $dna_start, $max - $min );
            }
            else {              #intergenic?

            }
        }

        my ( $codon, $codon_pos ) = get_codon( $t, $g, $mod_dna );
        my $nt = substr( $codon, $codon_pos - 1, 1 );
        die "ERROR - PAM position is not a G/C! "
          if lc($nt) ne 'g' and lc($nt) ne 'c';
        my $aa = $codons{ uc($codon) };
        my $gc = $nt eq 'g' ? 'c' : 'g';
        for my $m ( "A", uc($gc), "T" ) {
            my $mut_codon   = $codon;
            my $changed_cds = uc(substr( $mut_codon, $codon_pos - 1, 1, $m ));
            my $m_aa        = $codons{ uc($mut_codon) };
            if ( $m_aa eq $aa )
            { #hooray - we've got a synonymous change to the PAM site for our template
                my $g_m = $m;
                if ( $t->{strand} eq '-' ) {
                    $g_m = revcomp($g_m);
                }
                my $changed_base =
                  uc(substr( $mod_dna, $g - $dna_start, 1, uc($g_m) ));
                my $changed_c = get_coding_pos( $t, $g );
                my $distance_from_splice = get_exonic_distance_from_splice($t, $g);
                if ($distance_from_splice <= $splice_distance){
                    print "WARNING - exonic PAM site mutation is only $distance_from_splice away from a splice site!\n";
                }
                print
"Changing $t->{chrom}:$g from $changed_base to $g_m (c.$changed_c"
                  . $changed_cds
                  . ">$m) for synonymous PAM site mutation of repair template.\n";
                return substr( $mod_dna, $min - $dna_start, $max - $min );
            }
        }
    }
    return;
}

############################################
sub get_exonic_distance_from_splice{
    my ($t, $g_pos) = @_; 
    for ( my $i = 0 ; $i < @{ $t->{starts} } -1 ; $i++ ) {
        if ($g_pos >= $t->{starts}->[$i] and $g_pos <= $t->{ends}->[$i]){
            my $acc_distance = $g_pos - $t->{starts}->[$i] ;
            my $don_distance = $t->{ends}->[$i] - $g_pos + 1;
            return $acc_distance > $don_distance ? $don_distance : $acc_distance;
        }
    }
    die "Error getting $g_pos distance from splice site - is $g_pos exonic?\n";
}
############################################
sub get_intronic_position {
    my ( $t, $g_pos ) = @_;
    #returns coding coordinate and relative intronic coordinate
    for ( my $i = 0 ; $i < @{ $t->{starts} } -1 ; $i++ ) {
        if ( $g_pos >= $t->{ends}->[$i] && $g_pos < $t->{starts}->[$i+1] ) {
            my $intron_pos;
            my $exon_pos;
            my $donor_pos = $g_pos - $t->{ends}->[$i];
            my $acceptor_pos = $t->{starts}->[$i+1] - $g_pos + 1;
            if (abs($donor_pos) > abs($acceptor_pos)){
                $intron_pos =  $acceptor_pos;
                $exon_pos = $t->{starts}->[$i+1];
            }else{
                $intron_pos = $donor_pos;
                $exon_pos = $t->{ends}->[$i];
            }
            my $c_pos = get_coding_pos($t, $exon_pos); 
            
            if ( $t->{strand} eq '-' ) {
                $intron_pos *= -1;
            }
            return ($c_pos, $intron_pos); 
        }
    }
    die "Can't get intronic position for $g_pos for $t->{id} - is this coordinate really intronic?\n";
}
############################################
sub get_codon {
    my ( $t, $genomic_pos, $dna ) = @_;
    my $cdna = get_cdna( $t, $dna );
    my $c_pos = get_coding_pos( $t, $genomic_pos );

    #my $c = int( ($c_pos -1) /3) ; #0-based codon number
    my $p = $c_pos % 3;    #1-based position of $pos within codon
    $p ||= 3;
    my $codon_start = $c_pos - $p;    #0-based coordinate of beginning of codon
    my $codon = substr( $cdna, $codon_start, 3 );
    return ( $codon, $p );
}

############################################
sub is_intronic {
    my ( $t, $g_pos ) = @_;
    if ( $g_pos > $t->{ends}->[-1] ) {
        return 0;
    }
    if ( $g_pos < $t->{starts}->[0] ) {
        return 0;
    }
    return ( !is_exonic( $t, $g_pos ) );
}
############################################
sub is_exonic {
    my ( $t, $g_pos ) = @_;
    if ( $g_pos > $t->{ends}->[-1] ) {
        return 0;
    }
    if ( $g_pos < $t->{starts}->[0] ) {
        return 0;
    }
    for ( my $i = 0 ; $i < @{ $t->{starts} } ; $i++ ) {
        if ( $g_pos >= $t->{starts}->[$i] && $g_pos <= $t->{starts}->[$i] ) {
            return 1;
        }
    }
    return 0;
}
############################################
sub is_coding {
    my ( $t, $g_pos ) = @_;
    my @cds = get_cds($t);
    for ( my $i = 0 ; $i < @cds - 1 ; $i += 2 ) {
        if ( $g_pos >= $cds[$i] && $g_pos <= $cds[ $i + 1 ] ) {
            return 1;
        }
    }
    return 0;
}
############################################
sub get_utr_pos {
    my ($t, $g_pos) = @_; 
    my $u_pos;
    my $utr_length;  
    if ($g_pos < $t->{cdsStart}){
        for ( my $i = 0 ; $i < @{ $t->{starts} } ; $i++ ) {
            last if $t->{starts}->[$i] > $t->{cdsStart};#exon is in cds
            if ($t->{ends}->[$i] < $t->{cdsStart}){#whole exon is pre cds
                my $exon_length = $t->{ends}->[$i] - $t->{starts}->[$i];
                $utr_length += $exon_length;
                if ($g_pos > $t->{ends}->[$i]){
                    $u_pos += $exon_length;
                }elsif($g_pos >= $t->{starts}->[$i]){
                    $u_pos += $g_pos - $t->{starts}->[$i];
                }
            }else{#cds start is in exon
                $utr_length += $t->{cdsStart} - $t->{starts}->[$i];
                if ($g_pos >= $t->{starts}->[$i]){
                    $u_pos += $g_pos - $t->{starts}->[$i];
                }
                last; 
            }
        }
        $u_pos = $utr_length - $u_pos;
    }else{
        for ( my $i = 0 ; $i < @{ $t->{starts} } ; $i++ ) {
            next if $t->{ends}->[$i] < $t->{cdsEnd};#exon is before or within cds
            if ($t->{starts}->[$i] > $t->{cdsEnd}){#whole exon is post cds
                my $exon_length = $t->{ends}->[$i] - $t->{starts}->[$i];
                $utr_length += $exon_length;
                if ($g_pos > $t->{ends}->[$i]){
                    $u_pos += $exon_length;
                }elsif($g_pos >= $t->{starts}->[$i]){
                    $u_pos += $g_pos - $t->{starts}->[$i];
                }
            }else{#cds end is in exon
                $utr_length += $t->{ends}->[$i] - $t->{cdsEnd};
                if ($g_pos >= $t->{cdsEnd}){
                    $u_pos += $g_pos - $t->{cdsEnd};
                }
                last; 
            }
        }
    }
    if ( $t->{strand} eq '-' ) {
        $u_pos *= -1;
    }
    if ($u_pos > 0){
        $u_pos = "*$u_pos";
    }
    return $u_pos;
}
############################################
sub get_coding_pos {
    my ( $t, $g_pos ) = @_;
    return get_utr_pos($t, $g_pos) if not is_coding( $t, $g_pos );
    my @cds = get_cds($t);
    my $pos = 0;
    for ( my $i = 0 ; $i < @cds - 1 ; $i += 2 ) {
        last if ( $cds[$i] > $g_pos );
        if ( $cds[ $i + 1 ] < $g_pos ) {    #$g_pos is after this exon
            $pos += $cds[ $i + 1 ] - $cds[$i];
        }
        else {                              #g_pos is before end of exon
            if ( $g_pos > $cds[$i] ) {
                $pos += $g_pos - $cds[$i];
                last;
            }
        }
    }
    if ( $t->{strand} eq '-' ) {
        $pos = get_exons_length( \@cds ) - $pos + 1;
    }
    return $pos;
}
############################################
sub get_cdna {
    my ( $t, $dna ) = @_;
    my @cds  = get_cds($t);
    my $cdna = "";
    for ( my $i = 0 ; $i < @cds - 1 ; $i += 2 ) {
        $cdna .= substr(
            $dna,
            $cds[$i] - $t->{starts}->[0] + 1,
            $cds[ $i + 1 ] - $cds[$i],
        );
    }
    if ( $t->{strand} eq '-' ) {
        return revcomp($cdna);
    }
    return $cdna;

}

############################################
sub get_cds {
    my ($t) = @_;
    die "Transcript $t->{id} is non coding!\n"
      if $t->{cdsStart} == $t->{cdsEnd};
    my @cds = ();
    if ( exists $t->{cdsExonStarts} && @{ $t->{cdsExonStarts} } ) {
        for ( my $i = 0 ; $i < @{ $t->{cdsExonStarts} } ; $i++ ) {
            push @cds, $t->{cdsExonStarts}->[$i];
            push @cds, $t->{cdsExonEnds}->[$i];
        }
        return @cds;
    }
    push( @cds, $t->{cdsStart} );
    for ( my $i = 0 ; $i < @{ $t->{starts} } ; $i++ ) {
        next if ( $t->{ends}->[$i] < $t->{cdsStart} );
        push( @cds, $t->{cdsEnd} )
          if (  $t->{starts}->[$i] <= $t->{cdsStart}
            and $t->{ends}->[$i] >= $t->{cdsEnd} )
          ;    # if cds start and stop is in same exon
        push( @cds, $t->{starts}->[$i] )
          if (  $t->{starts}->[$i] > $t->{cdsStart}
            and $t->{starts}->[$i] < $t->{cdsEnd} )
          ; # get exonstart coord if greater than cds_start and less than cds end
        push( @cds, $t->{ends}->[$i] )
          if (  $t->{ends}->[$i] > $t->{cdsStart}
            and $t->{ends}->[$i] < $t->{cdsEnd} )
          ;  # get exonend coord if greater than cds start and less than cds end
        push( @cds, $t->{cdsEnd} )
          if (  $t->{ends}->[$i] >= $t->{cdsEnd}
            and $t->{starts}->[$i] > $t->{cdsStart} );
        last if ( $t->{ends}->[$i] >= $t->{cdsEnd} );
    }
    if ( @cds % 2 != 0 ) {
        die "Transcript $t->{id} contains an unresolved CDS.\n";
    }
    for ( my $i = 0 ; $i < @cds - 1 ; $i += 2 ) {

        #add to transcript hash as well for future reference
        push @{ $t->{cdsExonStarts} }, $cds[$i];
        push @{ $t->{cdsExonEnds} },   $cds[ $i + 1 ];
    }
    return @cds;
}

############################################
sub get_exons_length {
    my ($ex_array) = @_;
    my $length = 0;
    for ( my $i = 0 ; $i < @{$ex_array} - 1 ; $i += 2 ) {
        $length += ( $ex_array->[ $i + 1 ] - $ex_array->[$i] );
    }
    return $length;
}

############################################
sub get_cdna_genomic_coordinate {
    my ( $t, $p, $intron ) = @_;
    if ($intron) {
        $intron = -1 * $intron if $t->{strand} eq '-';
    }
    else {
        $intron = 0;
    }
    my $incr = 0;
    die "Transcript $t->{id} is non coding!\n"
      if $t->{cdsStart} == $t->{cdsEnd};
    my @cds           = get_cds($t);
    my $coding_length = get_exons_length( \@cds );
    die
"User coordinate $p is greater than the coding length ($coding_length) of $t->{id}.\n"
      if $coding_length < $p;
    my $cds_coord = 0;
    if ( $t->{strand} eq "+" ) {

        for ( my $i = 0 ; $i < @cds - 1 ; $i += 2 ) {
            $cds_coord += $cds[ $i + 1 ] - $cds[$i];
            if ( $p <= $cds_coord ) {
                my $diff = $cds_coord - $p;
                return $cds[ $i + 1 ] - $diff + $intron;
            }
        }
    }
    elsif ( $t->{strand} eq "-" ) {
        for ( my $i = $#cds ; $i > 0 ; $i -= 2 ) {
            $cds_coord += $cds[$i] - $cds[ $i - 1 ];
            if ( $p <= $cds_coord ) {
                my $diff = $cds_coord - $p;
                return $cds[ $i - 1 ] + $diff + 1 + $intron;
            }
        }
    }
    die "Error - genomic coordinate retrieval failed for $t->{id} c.$p\n";
}
############################################
sub get_gene_coordinates {
    my $t_id = shift;
    $t_id =~ s/\.\d+$//;    #remove version numbers if present
    my %t = ();             #details of exons and cds

    #CONNECT AND EXECUTE UCSC QUERY
    my $dbh = DBI->connect_cached( "dbi:mysql:$build:genome-mysql.cse.ucsc.edu",
        "genome", '' );
    my $command =
"SELECT name, chrom, cdsStart, cdsEnd, exonStarts, exonEnds, strand, name2 FROM $build.refGene WHERE name=?";
    my $sth = $dbh->prepare($command);
    $sth->execute($t_id);

    #PARSE RESULTS
    my $found = 0;
    while ( my @row = $sth->fetchrow_array ) {
        $found++;
        $t{id}       ||= $row[0];
        $t{symbol}   ||= $row[7];
        $t{strand}   ||= $row[6];
        $t{cdsStart} ||= $row[2];
        $t{cdsEnd}   ||= $row[3];
        $t{chrom}    ||= $row[1];
        @{ $t{starts} } = split( ",", $row[4] );
        @{ $t{ends} }   = split( ",", $row[5] );
    }
    die "No transcripts found for $t_id\n" if not $found;
    $t{dna} = get_transcript_dna( \%t );
    return %t;
}

############################################
sub get_transcript_dna {
    my $t   = shift;
    my $das = get "http://genome.ucsc.edu/cgi-bin/das/$build/dna?segment"
      . "=$t->{chrom}:$t->{starts}->[0],$t->{ends}->[-1]";
    my $xml  = new XML::Simple;
    my $data = $xml->XMLin($das);
    if ( not $data->{SEQUENCE}->{DNA}->{content} ) {
        die
"ERROR - No DNA found in $build for gene (chr$t->{chrom}:$t->{starts}->[0],$t->{ends}->[-1])!\n";
    }
    my $dna = lc( $data->{SEQUENCE}->{DNA}->{content} );
    $dna =~ s/\s//g;
    return $dna;
}
############################################
sub revcomp {
    my $seq     = shift;
    my $revcomp = reverse $seq;
    $revcomp =~ tr/acgtACGT/tgcaTGCA/;
    return $revcomp;
}
#############################################
sub search_sequence {
    my ( $dna, $search_sequence ) = @_;
    my @found = ();
    while ( $dna =~ /$search_sequence/ig ) {
        push( @found, pos($dna) + 1 );
    }
    return @found;
}

#############################################
sub parse_note {
    my $value = shift;
    my ( $score, $seq );
    if ( $value =~ /"score": "(\d+)%"/ ) {
        $score = $1;
    }
    if ( $value =~ /"sequence": "(\w+)"/ ) {
        $seq = $1;
    }
    return ( $score, $seq );
}

#############################################
sub parse_moiety {
    my $value = shift;
    my ( $strand, $id );
    if ( $value =~ /(\w+)-strand/ ) {
        $strand = $1;
    }
    if ( $value =~ /guide #(\d+)/ ) {
        $id = $1;
    }
    return ( $strand, $id );
}

#############################################
sub usage {
    my $msg = shift;
    if ($msg) {
        print "ERROR: $msg\n\n";
    }
    print <<EOT;
    Produces suitable guide sequence oligos and repair templates for homology directed repair using output from the Zhang lab CRISPR designing software http://crispr.mit.edu/.
    
    Usage: perl $0 -g <crispr genbank file> -t <refseq transcript id> -c <coding position> -m <mutation>  [options]
    
    Options:
    
    -t    --transcript      [RefSeq transcript ID of your target gene]
    -g    --genbank         [genbank file for guide RNAs produced by http://crispr.mit.edu/]
    -s    --score           [minimum guide RNA score to consider. Default = 50]
    -c    --coordinate      [cDNA position in transcript to mutate]
    -i    --intron_position [give this value in conjunction with -c to give intron/UTR position if your change is not coding]
    -m    --mutation        [mutation to introduce. Only SNVs supported.]
    -f    --flanks          [size of repair template flanks. Default = 60]
    -p    --forward_prefix  [sequence to prepend to forward guide. Default = CACC]
    -r    --reverse_prefix  [sequence to prepend to reverse guide. Default = AAAC]
    -e    --extra           [add this extra seq after the forward prefix to the forward guide and reverse complement to the 3' end of the reverse guide. Default = G]
    -d    --distance        [print a warning message if PAM site mutation is this many bp or less away from a splice site. Default is 8.]
    -b    --build           [genome build to use. Default = hg19]
    -h    --help            [show this help message and exit] 

Copyright 2015  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

EOT
    exit 1 if $msg;
    exit;
}


