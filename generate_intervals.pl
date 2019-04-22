#!/bin/env perl

use strict;
#use File::Basename;
#use lib dirname($0);
use lib '/home/ec2-user/tools/extra/tmp-work/github/bioinfo-tools';
use PM::Interval;

my @args = @ARGV;

if($#args < 0)
{
	die "
Usage: $0 <input-interval> [<bin-size>]

This programs reads an input file with genomic intervals and then 
output bins with the specified <bin-size>.

The input intervals should follow the format:
seq1<tab>start1<tab>end1
seq2<tab>start2<tab>end2
... ...

The coordinates are in 1-based right-closed normal mode.

The default <bin-size> is 300.

E.g.: $0 input.tsv 1000

";
}

my $file = shift;
my $binSize = shift;
$binSize ||= 300;

my %params = (size => $binSize, 'slice-interval' => 1);
#$params{'seqs'} = ['chr21','chrX'];

if(-f $file)
{
	$params{'file'} = $file;
}else
{
	$params{'build'} = $file;
}

my $handler = PM::Interval->new(%params);

my $counter= 0;
while(my $interval = $handler->next)
{
	printf(STDERR "%5d intervals have been output\n", $counter)
	if(++$counter % 100000 == 0);

	printf("%s\t%d\t%d\n",@$interval);
}

exit 0;

