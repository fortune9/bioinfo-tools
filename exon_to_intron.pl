#!/usr/bin/env perl

use strict;
use warnings;

# all exons in a transcript must appear in consecutive lines in the
# input file
my @exonSet; # store all exons in an mRNA
my ($rna, $exonNum);

my $counter = 0;

while(<>) {
    next if /^#/;
    chomp;
    my @fields = split "\t";
    if(++$counter % 50000 == 0 ) {
        warn "$counter lines are processed\n";
    }
    if(not $rna) { # first line
        ($rna, $exonNum) = get_rna_exon($fields[3]);
        $exonSet[0] = $rna;
        $exonSet[$exonNum] = \@fields;
        next;
    }
    # exons exist
    ($rna, $exonNum) = get_rna_exon($fields[3]);
    if($rna ne $exonSet[0]) { # new RNA
        my $intronsRef = extract_introns(\@exonSet);
        unless($intronsRef) {
            warn "No introns for $exonSet[0]\n";
        }
        output_introns($intronsRef);
        # re-initialize
        @exonSet = ();
        $exonSet[0] = $rna;
        $exonSet[$exonNum] = \@fields;
    } else { # add the exon to the RNA
        $exonSet[$exonNum] = \@fields;
    }
}

# don't forget the last one
if(@exonSet) {
    my $intronsRef = extract_introns(\@exonSet);
    output_introns($intronsRef);
}

warn "Job [$counter lines in total] is done\n";

exit 0;


sub get_rna_exon {
    my $str = shift;
    return split(':', $str);
}

sub extract_introns {
    my $exonsRef = shift;

    my $size = $#$exonsRef;
    if($size == 1) { # one exon, no intron
        return undef;
    }

    my @introns;
    my $seq = $exonsRef->[1]->[0];
    my $rna = $exonsRef->[0];
    for(my $i = 2; $i <= $size; $i++) {
       my $start = $exonsRef->[$i-1]->[2]; # end of previous 
       my $end = $exonsRef->[$i]->[1]; # start of current
       if($start > $end) { # minus strand
           $start = $exonsRef->[$i]->[2]; # end of current
           $end = $exonsRef->[$i-1]->[1]; # start of previous
       }
       # convert it to 0-based UCSC coordinates
       # $start -= 1; # start is the last base of previous exon,
       # giving the first base of the intron based on 0-based
       # coordinates.
       $end -= 1;
       my @intron = ($seq, $start, $end,$rna.':'.($i-1));
       my $tmpExon=$exonsRef->[$i-1];
       if($#$tmpExon >= 4) { # copy more columns from previous exon
           push @intron, $tmpExon->@[4..$#$tmpExon];
       }
       push @introns, \@intron;
    }
    return \@introns;
}

sub output_introns {
    my $intronsRef = shift;
    return unless $intronsRef; # No introns
    foreach (@$intronsRef) {
        print join("\t", @$_), "\n";
    }
}

