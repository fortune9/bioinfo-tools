#!/bin/bash

if [[ $# -lt 2 ]]; then
	cat <<EOF
Usage: $0 <bed-file1> <bed-file2> [<sorted?> [<output-prefix>]]

This program compares two bed-format region files and output three
files:

1. <output-prefix>.common.bed: the regions overlapped (see below),
which includes the input regions and their overlapped length.
2. <output-prefix>.1_specific.bed: the regions specific to <bed-file1>.
3. <output-prefix>.2_specific.bed: the regions specific to <bed-file2>.

The input files can be gzipped.
If the input files have been sorted by coordinates, feed the 3rd
argument with 'T', otherwise any other values to let the program to
sort.

The overlap between two regions is used to determine common regions
between the two files. Right now, if >= 50% of a region in <bed-file1>
is overlapped with any region in <bed-file2>, or vice versa, then the
region is considered common in botn input files.

The default value for output prefix is bed1_vs_bed2.\$\$, where \$\$
is the PID.

EOF
	exit 1;

fi

bedF1=$1;
bedF2=$2;
sorted=${3:-""}
outPre=${4:-bed1_vs_bed2.$$};
#echo  ">$sorted<", $outPre

if [[ "$sorted" != "T" ]]; then
	sorted=""
fi
# global variables
minFrac1=0.5;
minFrac2=0.5;
commFile=${outPre}.common.bed;
bed1Spec=${outPre}.1_specific.bed;
bed2Spec=${outPre}.2_specific.bed;

if [[ ! $(command -v bedtools) ]]; then
	echo "[ERROR] No command 'bedtools' found. Please install it"
	exit 2;
fi

# use a different variable for input files, because sorting is needed
# when necessary, without changing original files
if [[ ! $sorted ]]; then
	echo "[INFO] Sort input files in coordinates"
	inF1="inFile.1.$$.bed"
	inF2="inFile.2.$$.bed"
	less $bedF1 | sort -k1,1 -k2,2n >$inF1
	less $bedF2 | sort -k1,1 -k2,2n >$inF2
else
	inF1=$bedF1;
	inF2=$bedF2;
fi

echo "[INFO] Get overlapped regions into '$commFile'"
bedtools intersect -wo -a $inF1 -b $inF2 -f $minFrac1 -F $minFrac2 \
-e -sorted >$commFile;

echo "[INFO] Get $bedF1 specific regions into '$bed1Spec'"
bedtools intersect -v -a $inF1 -b $inF2 -f $minFrac1 -F $minFrac2 \
-e -sorted >$bed1Spec

echo "[INFO] Get $bedF2 specific regions into '$bed2Spec'"
bedtools intersect -v -b $inF1 -a $inF2 -f $minFrac1 -F $minFrac2 \
-e -sorted >$bed2Spec

echo "Job started at `date`"

# remove generated sorted files if applicable
if [[ ! $sorted ]]; then
	echo "Deleting temporary sorted files"
	rm $inF1;
	rm $inF2;
fi

echo "Job done at `date` !!!"

exit 0;

