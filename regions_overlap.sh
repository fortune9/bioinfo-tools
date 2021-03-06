#!/bin/bash

set -e

## functions

function cmp_num()
{
	if [[ $(echo "$1 $3 $2" | bc ) -gt 0 ]]; then 
		echo 1
	else
		echo ""
	fi
}

function summary_file
{
	echo "`gawk	'BEGIN{FS="\t";OFS="\t";i=0;s=0}{i++;s+=$3-$2;}END{print i,s;}'	<(less $1)`"
}

function warn
{
	echo $* >&2
}

if [[ $# -lt 2 ]]; then
	cat <<EOF
Usage: $0 <bed-file1> <bed-file2> [<frac1> <frac2> <sorted?> [<output-prefix>]]

This program compares two bed-format region files and output these
files:

1. <output-prefix>.common.bed: the regions overlapped (see below),
which includes the input regions and their overlapped length.
2. <output-prefix>.1_specific.bed: the regions specific to <bed-file1>.
3. <output-prefix>.2_specific.bed: the regions specific to <bed-file2>.
4. <output-prefix>.summary.tsv: gives a summary of results.

The input files can be gzipped.

The common regions are the ones overlap between <bed-file1> and
<bed-file2>. Here 'overlap' means the following 2 conditions are
all satisfied:

1. The overlapped region covers <frac1> fraction or more of the
involved region from <bed-file1>.

2. The overlapped region covers <frac2> fraction or more of the
involved region from <bed-file2>.

Otherwise the region is classified as file-specific.

One can specify a value between 0 and 1 (inclusive) for the arguments
<frac1> and <frac2> to define the overlap requirement. When the value
is 0, it means the overlap requirement is 1bp, which is also the
default.

If the input files have been sorted by coordinates, feed the 5th
argument with 'T', otherwise any other values to let the program to
sort.


The default value for output prefix is bed1_vs_bed2.\$\$, where \$\$
is the PID.

EOF
	exit 1;

fi

bedF1=$1;
bedF2=$2;
minFrac1=${3:-1e-10}
minFrac2=${4:-1e-10}
sorted=${5:-""}
outPre=${6:-bed1_vs_bed2.$$};
#echo  ">$sorted<", $outPre

if [[ "$sorted" != "T" ]]; then
	sorted=""
fi

if [[ $(cmp_num "$minFrac1" 0 "==") ]]; then
	minFrac1=1e-10;
fi

if [[ $(cmp_num "$minFrac2" -0 "==") ]]; then
	minFrac2=1e-10;
fi

#echo "I got frac1:'$minFrac1', frac2: '$minFrac2'"; exit 0;

# global variables
commFile=${outPre}.common.bed;
bed1Spec=${outPre}.1_specific.bed;
bed2Spec=${outPre}.2_specific.bed;
sumFile=${outPre}.summary.tsv;

if [[ ! $(command -v bedtools) ]]; then
	warn "[ERROR] No command 'bedtools' found. Please install it"
	exit 2;
fi

# use a different variable for input files, because sorting is needed
# when necessary, without changing original files
if [[ ! $sorted ]]; then
	warn "[INFO] Sort input files in coordinates"
	inF1="inFile.1.$$.bed"
	inF2="inFile.2.$$.bed"
	less $bedF1 | sort -k1,1 -k2,2n >$inF1
	less $bedF2 | sort -k1,1 -k2,2n >$inF2
else
	inF1=$bedF1;
	inF2=$bedF2;
fi

warn "Job started at `date`"

warn "[INFO] Get overlapped regions into '$commFile'"
bedtools intersect -wo -a $inF1 -b $inF2 -f $minFrac1 -F $minFrac2 \
-e -sorted >$commFile;

warn "[INFO] Get $bedF1 specific regions into '$bed1Spec'"
bedtools intersect -v -a $inF1 -b $inF2 -f $minFrac1 -F $minFrac2 \
-e -sorted >$bed1Spec

warn "[INFO] Get $bedF2 specific regions into '$bed2Spec'"
bedtools intersect -v -b $inF1 -a $inF2 -f $minFrac1 -F $minFrac2 \
-e -sorted >$bed2Spec

warn "[INFO] Summarize the results in '$sumFile'"

total1Sum=$(summary_file $inF1)
total2Sum=$(summary_file $inF2)
spec1Sum=$(summary_file $bed1Spec )
spec2Sum=$(summary_file $bed2Spec )
file1NF=$(less $inF1 | head -1 | gawk '{print NF}')
tmp1F=tmp1.$$
tmp2F=tmp2.$$
less $commFile | cut -f 1-$file1NF | uniq >$tmp1F
file1NF=$(( file1NF + 1 ))
less $commFile | cut -f $file1NF- | cut -f 1-4 | sort | uniq >$tmp2F
comm1Sum=$(summary_file $tmp1F)
comm2Sum=$(summary_file $tmp2F)

cat >$sumFile <<EOF
Type	Count	Size
1_all	$total1Sum
2_all	$total2Sum
1_comm	$comm1Sum
2_comm	$comm2Sum
1_spec	$spec1Sum
2_spec	$spec2Sum
EOF

rm $tmp1F $tmp2F;

# remove generated sorted files if applicable
if [[ ! $sorted ]]; then
	warn "Deleting temporary sorted files"
	rm $inF1;
	rm $inF2;
fi

warn "Job done at `date` !!!"

exit 0;

