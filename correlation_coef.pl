#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $pg = $0;
$pg =~ s!.*/!!;

my $sep = "\t";

my $cols;
my $noZero=0;
my $pseudoCnt;

GetOptions(
	"col:s"	=> \$cols,
	"sep:s"	=> \$sep,
	"no-zero!" => \$noZero,
	"log2:f" => \$pseudoCnt
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
my $nForLog2 = 0; # counts for log2-value pairs
my $xSum = 0;
my $xSqSum = 0;
my $ySum = 0;
my $ySqSum = 0;
my $xySum = 0;
my $xLog2Sum = 0;
my $xSqLog2Sum = 0;
my $yLog2Sum = 0;
my $ySqLog2Sum = 0;
my $xyLog2Sum = 0;

my $counter = 0;

if($noZero)
{
	warn "[$pg] Data pairs being zeros will be excluded\n";
}

while(my $values = next_value_pair())
{
	my ($x, $y) = @$values;
	if($noZero and $x==0 and $y==0) { next; }
	$xSum += $x;
	$xSqSum += $x**2;
	$ySum += $y;
	$ySqSum += $y**2;
	$xySum += $x*$y;
	$n++;
	if(defined $pseudoCnt)
	{
		$x+=$pseudoCnt;
		$y+=$pseudoCnt;
		if($x <= 0 or $y<=0) {next;}
		$x=log2($x); $y=log2($y);
		$xLog2Sum += $x;
		$xSqLog2Sum += $x**2;
		$yLog2Sum += $y;
		$ySqLog2Sum += $y**2;
		$xyLog2Sum += $x*$y;
		$nForLog2++;
	}
	warn "[$pg] $counter data points have been processed\n"
	if(++$counter % 100000 == 0);
}

my $xSig=sqrt($n*$xSqSum - $xSum**2);
my $ySig=sqrt($n*$ySqSum - $ySum**2);
my $r = "NA";
my $s="%s";

printf "#%s\t%s\t%s\t%s\t%s\t%s\n", "sample_size", "r", "mean_x", "sigma_x", "mean_y", "sigma_y";
if($xSig != 0 and $ySig !=0)
{
	$r = ($n*$xySum - $xSum*$ySum)/($xSig*$ySig);
	$s="%.6g";
}else
{
	warn "The variances for input variables are 0; can't compute\n";
	#	printf "%d\tNA\t%.6g\t%.6g\t%.6g\t%.6g\n", $n, $xSum/$n,$xSig/$n,$ySum/$n,$ySig/$n;
}
printf "%d\t$s\t%.6g\t%.6g\t%.6g\t%.6g\n", 
	$n, $r, $xSum/$n,sample_var($n,$xSig),$ySum/$n,sample_var($n,$ySig);

if($nForLog2 > 0)
{
	printf "#%s\t%s\t%s\t%s\t%s\t%s\t(For log2 data)\n", 
	"sample_size", "r", "mean_x", "sigma_x", "mean_y", "sigma_y";
	$xSig=sqrt($nForLog2*$xSqLog2Sum - $xLog2Sum**2);
	$ySig=sqrt($nForLog2*$ySqLog2Sum - $yLog2Sum**2);
	$r='NA';
	$s="%s";
	if($xSig != 0 and $ySig !=0)
	{
		$r = ($nForLog2*$xyLog2Sum - $xLog2Sum*$yLog2Sum)/($xSig*$ySig);
		$s="%.6g";
		#	printf "%d\t%.6g\n", $nForLog2, $r;
	}else
	{
		warn "The variances for log2-transformed variables are 0; can't compute\n";
		#printf "%d\tNA\n", $nForLog2;
	}
	
	printf "%d\t$s\t%.6g\t%.6g\t%.6g\t%.6g\n", 
		$nForLog2, $r, $xLog2Sum/$nForLog2,sample_var($nForLog2,$xSig),
		$yLog2Sum/$nForLog2,sample_var($nForLog2, $ySig);

}

close $fh1;
if(defined $inFile2) { close $fh2; }

warn "[$pg] Job is done at ".localtime()."\n";

exit 0;

sub sample_var
{
	my ($n, $nSigma) = @_;
	return $nSigma/sqrt($n*($n-1)); # need correct sample size
}

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

sub log2
{
	my $x = shift;
	return log($x)/log(2);
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

--log2: <num>, if this option is provided, then correlation for
log2-transformed data will also be computed, and the provided value
for this option will be used as a pseudocount to add to all values
before transformation. []

--sep:  <string>, the field separator for input files. [tab]

--no-zero: a switch option. When provided, the data point pairs
whose values are all zeros are excluded in calculation. This is
useful when computing read depths correlation between two samples.

Example uses:
# value in 2 files
$0 --col 1,3 x.txt y.txt
# value in 1 file
$0 --col 1,3 all.txt
# value in 1 file and exclude rows with all zeros
$0 --no-zero --col 1,3 all.txt
# also compute correlation for log2-transformed data with
# pseudocount=1
$0 --no-zero --col 1,3 --log2 1  all.txt

EOF

}

