#!/bin/bash

function usage()
{
	cat <<EOF
Usage: $0 <input-file> [<field-sep>]

This program reads a file and transpose it, columns -> rows.

The option <field-sep> gives the field separator for input file,
default is tab.

E.g.: $0 columns.tsv "\t"

EOF

}

function warn()
{
	echo $* >&2
}

if [[ $# -lt 1 ]]; then
	usage;
	exit 1;
fi

f=$1;
sep=$2;

if [[ "$sep" == "" ]]; then
	sep="\t"
fi

nf=`head -1 $f | gawk -v s=$sep 'BEGIN{FS=s}{print NF}'`

for i in `seq 1 $nf`;
do
	if [[ "$sep" == "\t" ]]; then
		fs=""
	else
		fs="-d $sep"
	fi
	newRow=$(cut -f $i $fs $f | tr '\n' "$sep")
	#newRow=${newRow%$sep}
	echo "$newRow" | sed -e "s/$sep$//"
done

warn "Job is done"

exit 0;

