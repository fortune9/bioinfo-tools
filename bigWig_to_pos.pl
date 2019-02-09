#!/usr/bin/env perl
use strict;

my $file = shift or die "Usage: $0 <bigwig file>:$!";

open(o, "bigWigToWig $file stdout | ") or die $!;

my $chrom="";
my $span=0;

my $counter = 0;
my $mode="vs"; # default is variable step
my $firstLine = <o>;
if($firstLine =~ /^#/) { # bedgraph mode
	$mode='bg';
}else{
	output_vs_line($firstLine); # a vs line
}
warn "#working mode is $mode for $file\n";

while(<o>) {
	chomp;
	if($mode eq 'bg') # for data from IHEC
	{
		next if /^#/;
		my @fields = split "\t";
		die "The region length in '$_' is not 2:$!" if $fields[2]-$fields[1] != 2;
		$fields[1]+=1; # make it 1-based format
		print join("\t", @fields[0,1,3]),"\n";
		$counter++;
		next;
	}
	# otherwise the default variableStep mode
	output_vs_line($_);
	$counter++;
}

close o;

warn "# totally $counter lines in $file are processed\n";

exit 0;

sub output_vs_line
{
	my $line = shift;
	if($line=~/^variableStep\s+chrom=(\S+)\s+span=(\d+)/i)
	{
		$chrom=$1; $span=$2;
		print STDERR join("\t", $chrom, $span),"\n";
	}else {
		my ($pos, $d) = split /\s+/;
		printf("%s\t%d\t%d\n", $chrom, $pos, $d); 
	} 
}
