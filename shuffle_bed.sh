#!/bin/bash

set -e

# load common functions
bashFuncUrl='https://raw.githubusercontent.com/fortune9/programs/master/bash_functions.sh'
bashFuncFile="tmp.$$.bash_functions.sh"

wget --quiet -O $bashFuncFile "$bashFuncUrl"

if [[ ! -f $bashFuncFile ]]; then
	echo "Can't download 'bash_functions.sh' from github"
	exit 5;
fi

source $bashFuncFile


function usage()
{
	cat <<EOF
Usage: $0 [options] <genome-file> <bed-file> [<extra args directly go
to bedtools>]

This program shuffle the genomic intervals in the <bed-file> along the
genome given by the <genome-file>. It can generate multiple batches.

The <bed-file> can be in gzipped format and can contain more than 3
columns;a unique id at the 4th column is recommended as this can be
used to track which interval is sampled in the output file.

<genome-file> is a tab-delimited file with chromosome name and length
at 1st and 2nd column, respectively.

All input bed files, including those for '--include' and '--exclude'
can be gzipped.

Options (default values are in []):

-h|--help   print out this usage information.

--include <path>   a bed file. If given, then the generated bed
intervals will be only from these bed intervals. []

--exclude <path>  a bed file. If given, then generated bed intervals
will not overlap any of the intervals in this file. []

--times <int>    an integer specifying the number of batches for
shuffling the input bed file. [1]

--out <path>   a file to store output. [shuffled.<bed-file>]

*Note*: when the option --include is specified, then the direct option
'-chrom' will not be respected by bedtools (a bug), so the use of
'--exclude' is preferred.

Example usages:

# shuffle a bed file 10 times
$0 --times 10 chr_size.tsv in.bed.gz

# shuffle the intervals within the same chromesome
$0 --times 10 chr_size.tsv in.bed.gz -chrom

# write to outfile 'out.bed'
$0 --times 10 --out out.bed chr_size.tsv in.bed.gz

# some tests
$0 --times 2 --out out1.bed chr_size.tsv in.bed.gz -chrom -seed 100
$0 --times 2 --out out2.bed --include incl.bed chr_size.tsv in.bed.gz -chrom -seed 100
$0 --times 2 --out out3.bed --exclude excl.bed chr_size.tsv in.bed.gz -chrom -seed 100

EOF
}

depends=(bedtools)

for e in ${depends[@]}
do
	if [[ ! $(check_exe $e) ]]; then
		echo "Command '$e' does not exist"
		exit 1;
	fi
done

if [[ $# -lt 2 ]]; then
	usage;
	exit 2;
fi

# read into parameters
inclBed="";
exclBed="";
times=1;
outFile="";
posArgs=();

while [[ $# -gt 0 ]];
do
	k=$1; shift;
	case $k in
		-h|--help)
			usage;
			exit 3;
			;;
		--include)
			inclBed=$1;
			shift;
			;;
		--exclude)
			exclBed=$1;
			shift;
			;;
		--times)
			times=$1;
			shift;
			;;
		--out)
			outFile=$1;
			shift;
			;;
		*) # all other inputs
			posArgs+=("$k")
			;;
	esac
done

# start work
if [[ ${#posArgs[@]} -lt 2 ]]; then
	echo "Both <genome-file> and <bed-file> inputs are needed"
	exit 2;
fi

set -- ${posArgs[@]}

chrSizeFile=$1; shift;
bedFile=$1; shift;
directArgs="$*"

#echo "$chrSizeFile : $bedFile : $directArgs"

if [[ $outFile == "" ]]; then
	outFile=shuffled.${bedFile%.gz}
fi

touch $outFile;
echo "" >$outFile;
exclStr="";
inclStr="";
maxFrac=1e-9;

if [[ $inclBed ]]; then inclStr="-incl $inclBed"; fi
if [[ $exclBed ]]; then exclStr="-excl $exclBed"; fi

for i in `seq $times`
do
	echo "Shuffling $bedFile: $i th time"
#	cat <<EOF
	bedtools shuffle -i $bedFile -g $chrSizeFile $inclStr $exclStr \
		-f $maxFrac -noOverlapping $directArgs \
		| gawk -v i=$i 'BEGIN{OFS="\t"}{print $0, "batch." i}' \
	>>$outFile
#EOF
done

echo "Job done and output is in '$outFile'"

# clean up
rm $bashFuncFile

exit 0;

