#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $type;
GetOptions(
    "type=s"  => \$type
);

my @fields;
my $counter = 0;
while(<>)
{
    next if /^#/;
    @fields = split "\t";
    if( ++$counter % 50000 == 0) {
        warn "$counter lines have been processed\n";
    }
    if($fields[2] ne $type) { next; }
    my %attrs = parse_attrs($fields[8]);
    if($type eq 'exon') {
        # seq, start, end, exon_id, score, strand, gene_id
        my $id = sprintf("%s:%d", $attrs{'transcript_id'}, $attrs{"exon_number"});
        print(join("\t", $fields[0], $fields[3]-1, $fields[4],
            $id, $fields[5], $fields[6], $attrs{'gene_id'}), "\n");
    } else {
        # seq, start, end, gene, score, strand
        print(join("\t", $fields[0], $fields[3]-1, $fields[4],
            $attrs{'gene_id'}, $fields[5], $fields[6]), "\n");
    }
}

warn "Job [$counter lines] is done\n";

exit 0;

sub parse_attrs {
    my $str = shift;
    chomp $str;
    # use quotes to separate the fields, two types: one followed by';'
    # and the other not
    my @data = split /\s*(?:"(?!;)|";)\s*/, $str;
    if(@data % 2 != 0) {
        die "Unknown data string: '$str', $!";
    }
    my %res = @data;
    return %res;
}

