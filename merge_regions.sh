#!/bin/bash

function usage
{
	cat <<EOF
Usage: $0 <bed-file1> <bed-file2> [<out-file>]

This program compares the regions from the two bed files and outputs
a file containing all the intervals overlapped and specific to each of
the inputs.

If the 3rd argument <out-file> is omitted, the result will be on the
standard output, normally, the screen.

EOF
	
}

depends=(bedtools)

for e in ${depends[@]}
do
	if [[ ! $(command -v $e) ]]; then
		echo "Command '$e' is not found"
		exit 2;
	fi
done

if [[ $# -lt 2 ]]; then
	usage;
	exit 1;
fi

f1=$1;
f2=$2;
o=$3;

tmpComm=tmp.$$.comm.bed
tmpAspec=tmp.$$.A.bed
tmpBspec=tmp.$$.B.bed

# get the common and specific intervals
bedtools intersect -a $f1 -b $f2 | gawk \
'BEGIN{OFS="\t";i=0}{i++; if(NF>3) {$4="comm.i" i "\t" $4;} else {$4="comm.i" i;} print $0}'	>$tmpComm
bedtools subtract  -a $f1 -b $f2 | gawk \
'BEGIN{OFS="\t"}{i++; if(NF>3) { $4= "A.spec.i" i "\t" $4; } else {$4="A.spec.i" i; } print $0}'	>$tmpAspec
bedtools subtract  -b $f1 -a $f2 | gawk \
'BEGIN{OFS="\t"}{i++; if(NF>3) { $4="B.spec.i" i "\t" $4; } else { $4="B.spec.i" i;} print $0}'	>$tmpBspec

if [[ "$o" ]]; then
	cat $tmpComm $tmpAspec $tmpBspec | sort -k1,1 -k2,2n > $o
else
	cat $tmpComm $tmpAspec $tmpBspec | sort -k1,1 -k2,2n
fi

rm $tmpComm $tmpAspec $tmpBspec

echo "Mergine is done" >&2

exit 0;

