#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $sep = "\t";

my $cols;

GetOptions(
	"col:s"	=> \$cols,
	"sep:s"	=> \$sep
);

my $inFile1 = shift;
my $inFile2 = shift;

unless($inFile1) { usage(); exit 1; }

my $fh1 = _open_file($inFile1);
my $fh2;
my $singleFile=0;

if(defined $inFile2)
{
	$fh2 = _open_file($inFile2);
	unless(defined $cols) { $cols="1,1"; }
}else
{
	$fh2=$fh1;
	warn "x and y values will be obtained from the same file\n";
	$singleFile=1;
	unless(defined $cols) { $cols="1,2"; }
}

my @cols=split ',', $cols;
@cols = map { $_-1 } @cols;

# now compute the sums of x and x^2
my $n = 0;
my $xSum = 0;
my $xSqSum = 0;
my $ySum = 0;
my $ySqSum = 0;
my $xySum = 0;

my $counter = 0;

while(my $values = next_value_pair())
{
	my ($x, $y) = @$values;
	$xSum += $x;
	$xSqSum += $x**2;
	$ySum += $y;
	$ySqSum += $y**2;
	$xySum += $x*$y;
	$n++;
	warn "$counter data points have been processed\n"
	if(++$counter % 100000 == 0);
}

my $xSig=sqrt($n*$xSqSum - $xSum**2);
my $ySig=sqrt($n*$ySqSum - $ySum**2);
my $r = "NA";

printf "#%s\t%s\n", "sample_size", "r";
if($xSig != 0 and $ySig !=0)
{
	$r = ($n*$xySum - $xSum*$ySum)/($xSig*$ySig);
	printf "%d\t%.6g\n", $n, $r;
}else
{
	warn "The variances for input variables are 0; can't compute\n";
	printf "%d\tNA\n", $n;
}


close $fh1;
if(defined $inFile2) { close $fh2; }

warn "Job is done at ".localtime()."\n";

exit 0;

sub _open_file
{
	my $f = shift;
	my $fh;
	if($f =~ /\.gz$/i)
	{
		open($fh, "zcat $f |") or die "Can't open gzipped $f:$!";
	}else
	{
		open($fh, "< $f") or die "Can't open $f:$!";
	}
	return $fh;
}

sub next_value_pair
{
	return undef if(eof($fh1) or eof($fh2));
	my $line1 = <$fh1>;
	chomp $line1;
	my @fields = split($sep, $line1);
	my $x=$fields[$cols[0]];
	my $y;
	if($singleFile)
	{
		$y=$fields[$cols[1]];
	}else
	{
		my $line2 = <$fh2>;
		chomp $line2;
		$y=(split($sep, $line2))[$cols[1]];
	}
	return undef unless(defined($x) and defined($y));
	return [$x,$y];
}

sub usage
{
	print <<EOF;
$0 [options] <value-file-1> [<value-file-2>]

This program computes Pearson's correlation coefficient 
using the formula based on the mean of x and sum of x^2.
Thus this is memory efficient. If the sample size is small,
using R's cor is recommended.

Here <value-file-1> and <value-file-2> specify the x and
y values for computing the correlation efficient. If only
<value-file-1> is provided, then both x and y values are
expected from this input file.

Input files can be gzipped; if so, please keep the extension '.gz'
in filenames.

Options (default values are in []):

--col:  <int1,int2>, the columns from which the x and y values
are present in the input file(s). The column number for a file
starts from 1. Default: [1,1] for two input files, [1,2] for
one input file.

--sep:  <string>, the field separator for input files. [tab]

Example uses:
# value in 2 files
$0 --col 1,3 x.txt y.txt
# value in 1 file
$0 --col 1,3 all.txt

EOF

}

