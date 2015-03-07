#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use DBI;
use LWP::Simple;
use XML::Simple;

my $pos;    #position in seq to mutate
my $int_pos = 0;    #intron/UTR position if not coding
my $transcript;     #target transcript RefSeq id
my $flanks = 30;
my $build  = 'hg19';
my $help;
GetOptions(
    "transcript=s"      => \$transcript,
    "coordinate=i"      => \$pos,
    "intron_position=i" => \$int_pos,
    "flanks=i"          => \$flanks,
    "build=s"           => \$build,
    "help"              => \$help,
) or usage("Error retriveing user options.");

usage()                                       if $help;
usage("--transcript argument is required.\n") if not $transcript;
usage("--coordinate argument is required.\n") if not $pos;

my %transcript = get_gene_coordinates($transcript);
my $pos_coord = get_cdna_genomic_coordinate( \%transcript, $pos, $int_pos );
print
"$transcript{id}: c.$pos matches genomic coordinate $transcript{chrom}:$pos_coord\n";
my $dna =
  get_dna( $transcript{chrom}, $pos_coord - $flanks, $pos_coord + $flanks );
my $format_dna = lc( substr( $dna, 0, $flanks ) );
$format_dna .= uc( substr( $dna, $flanks, 1 ) );
$format_dna .= lc( substr( $dna, $flanks + 1 ) );
print ">$transcript{id}|$transcript{symbol}|c.$pos";

if ($int_pos) {
    if ( $int_pos < 0 ) {
        print "$int_pos";
    }
    else {
        print "+$int_pos";
    }
}
print "\n$format_dna\n";

############################################
sub get_dna {
    my ( $chrom, $start, $end ) = @_;
    my $das = get "http://genome.ucsc.edu/cgi-bin/das/$build/dna?segment"
      . "=$chrom:$start,$end";
    my $xml  = new XML::Simple;
    my $data = $xml->XMLin($das);
    if ( not $data->{SEQUENCE}->{DNA}->{content} ) {
        die "ERROR - No DNA found for build $build ($chrom:$start,$end)!\n";
    }
    my $dna = lc( $data->{SEQUENCE}->{DNA}->{content} );
    $dna =~ s/\s//g;
    return $dna;
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
    my $dbh = DBI->connect_cached( "dbi:mysql:hg19:genome-mysql.cse.ucsc.edu",
        "genome", '' );
    my $command =
"SELECT name, chrom, cdsStart, cdsEnd, exonStarts, exonEnds, strand, name2 FROM hg19.refGene WHERE name=?";
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

    #$t{dna} = get_transcript_dna(\%t);
    return %t;
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

#############################################
sub usage {
    my $msg = shift;
    if ($msg) {
        print "ERROR: $msg\n\n";
    }
    print <<EOT;
    Usage: $0 -t <refseq transcript id> -c <coding position> [options]
    
    Options:
    
    -t    --transcript      [RefSeq transcript ID of your target gene]
    -c    --coordinate      [cDNA position in transcript to mutate]
    -i    --intron_position [give this value in conjunction with -c to give intron/UTR position if your change is not coding]
    -f    --flanks          [amount of bp flanks to retrieve either side of cDNA coordinate]
    -h    --help            [show this help message and exit] 

EOT
    exit 1 if $msg;
    exit;
}

