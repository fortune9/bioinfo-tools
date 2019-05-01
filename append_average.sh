#!/bin/bash

set -e

function msg_info
{
	echo "[INFO] $1" >&2
}

function msg_warn
{
	echo "[WARN] $1" >&2
}

function usage
{
	cat <<EOF
Usage: $0 <bed-file> <bedgraph-file> [<out-file>]

This program reads a bed file and calculates the average on the values
recorded on the 4th column of the <bedgraph-file>. Both input files
should be sorted in coordinates beforehand, and can be compressed.

The output is written to <out-file>, if provided; otherwise to
standard screen.

E.g.: $0 in.bed in.bedgraph.gz out.bed

EOF

}

if [[ $# -lt 2 ]]; then
	usage;
	exit 1;
fi

depends=(weighted_avg.py bedtools)

for e in "${depends[@]}"
do
	if [[ ! $( command -v $e ) ]]; then
		msg_warn "Command '$e' can not be found"
		exit 4;
	fi
done

bedF=$1;
bgF=$2;
outF=$3;

nf1=$( less $bedF | head -5 | gawk 'BEGIN{FS="\t"}NR==3{print NF}' )
nf2=$( less $bgF | head -5 | gawk 'BEGIN{FS="\t"}NR==3{print NF}' )

if [[ $nf2 -lt 4 ]]; then
	msg_warn "The bedgraph file '$bgF' has fewer than 4 fields"
	exit 3;
fi

if [[ $nf2 -ne 4 ]]; then
	msg_warn "The input bedgraph file '$bgF' has more than 4 fields"
fi

# the -wao garantee the regions w/o overlaps have a value 0 in output.
# we also add one more column in $bedF to give an unique id for each
# line
bedtools intersect -wao -a <(less $bedF | gawk \
'BEGIN{OFS="\t";i=0}{i++;l="Line." i; print $0, l}') -b $bgF -sorted >tmp1.$$

gf=$(( nf1 + 1 )); # group field: the uniq line id
vf=$(( nf1 + 5 )); # assume the 4th field for values
wf=$(( nf1 + nf2 + 2 ));

weighted_avg.py  -g $gf -v $vf -w $wf -o tmp2.$$ <(less tmp1.$$ | gawk \
-v v=$vf -v w=$wf 'BEGIN{OFS="\t"}{l=$3-$2; $w/=l; if($v==".") {$v=0;} print $0}' )

if [[ "$outF" ]]; then
	paste <(less -f $bedF) <(cut -f 2 tmp2.$$) > $outF
else
	paste <(less -f $bedF) <(cut -f 2 tmp2.$$) # to standard screen
fi

rm tmp1.$$ tmp2.$$

msg_info "Job is done"

exit 0;

