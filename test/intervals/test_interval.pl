#!/bin/env perl

use strict;
use File::Basename;
use lib dirname($0)."/../..";
use PM::Interval;

my $binSize = 30000;

my $file = shift  or die "$0 <file or species>\n";

my %params = (size => $binSize, 'slice-interval' => 1);
$params{'seqs'} = ['chr21','chrX'];

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
	printf("%5d:\t%s\t%d\t%d\n", ++$counter, @$interval);
}

exit 0;

