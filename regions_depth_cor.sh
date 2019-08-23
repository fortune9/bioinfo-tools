#!/bin/bash

set -e

tmpFiles=();

# load the common bash file
bashFuncUrl='https://raw.githubusercontent.com/fortune9/programs/master/bash_functions.sh'

bashFuncFile="tmp.$$.bash_functions.sh"
wget --quiet -O $bashFuncFile "$bashFuncUrl"
#cp /home/ubuntu/tools/extra/tmp-work/github/programs/bash_functions.sh $bashFuncFile

if [[ ! -f $bashFuncFile ]]; then
        echo "Can't download 'bash_functions.sh' from github"
        exit 5;
fi

source $bashFuncFile

tmpFiles+=("$bashFuncFile")

pg=${BASH_SOURCE[0]}
pgDir=`dirname $pg`

# add the program folder and pwd to PATH
PATH=.:$pgDir:$PATH;

# functions
## add average signal for a bed file and keep the input order
function add_signal
{
	bedFile=$1;
	sigFile=$2;
	nf=$(file_nf $bedFile)
	if [[ ! $(is_int $nf) ]]; then
		error "Unkown file format in '$bedFile': $nf"
		exit 3
	fi
	# sort this file
	tmpFile=$(rand_str).bed
	# add the line id and sort
	gawk 'BEGIN{FS="\t";OFS="\t";i=0}{i++; print $0, i}' $bedFile | \
		sort -k1,1 -k2,2n >$tmpFile
	# here assume the output keeps the order of the input regions by
	# bedtools.
	bedtools coverage -a $tmpFile -b $sigFile -hist -sorted \
	| gawk '$1 != "all"' >$tmpFile.cov
	gCol=$((nf + 1)) # line id
	vCol=$((gCol + 1))
	wCol=$((gCol + 4))
	python3 $(command -v weighted_avg.py) -g $gCol -v $vCol -w $wCol \
		$tmpFile.cov >$tmpFile.sig
	# combine and remove the lind id column
	paste <(cut -f 1-$nf $tmpFile) $tmpFile.sig | sort \
	-k$gCol,${gCol}n | gawk \
	'BEGIN{OFS="\t"}{pv=NF-1;$pv=$NF;NF=pv;print $0}' >$bedFile
	# clean up
	rm $tmpFile*
	echo $bedFile
}


function usage()
{
	cat <<EOF
Usage: $0 [options] <read-file1> <read-file2> <bed-file1> [<bed-file2>]

*Note*: all input files are allowed to be gzipped.

This program computes the Pearson correlation for the read depths
as well as for log2-transformed depths represented in the two 
mapped-reads files <read-file1> and <read-file2>. 

When only <bed-file1> is provided, it calculates read depths for the
regions in <bed-file1> using the <read-file1>, and then read depths
for the same regions using the <read-file2>, and finally computes the
correlation between the two sets of read depths.

When both <bed-file1> and <bed-file2> are provided, a set of read
depths is computed for regions in <bed-file1> using <read-file1> and
the other set of read depths in <bed-file2> using <read-file2>, and
then correlation is computed between the two sets of read depths. In
the scenario, the number of regions in the two bed files should be the
same and in comparable order.

<read-file1> and <read-file2> can be in bed/bam/gff format. If in
bed/gff format, they need be pre-sorted in coordinates with 'sort
-k1,1 -k2,2n'.

<bed-file1> and <bed-file2> needn't be pre-sorted.


--keep-avg: a switch option, when provided, the files containing
computed read depths for each region are kept in a file.

--no-zero: a switch option. If provided, then regions (or pairs)
with all samples being zeros are excluded.

--pseudo-count <num>: a number added to original average depth values
to compute log2-transformed data. This is necessary to save regions
with depth being zeros. [0]

Example uses:

# get the correlation between sample1 and sample2 for the regions in
# common.bed
$0 sample1.bam sample2.bam common.bed

# same as above, but keep the middle file with computed average depths
$0 --keep-avg sample1.bam sample2.bam common.bed

# exclude regions where read depths in both samples are 0
$0 --no-zero sample1.bam sample2.bam common.bed

EOF

rm $bashFuncFile

}

depends=(python3 weighted_avg.py correlation_coef.pl bam2bed.sh)

for i in ${depends[@]}
do
	if [[ ! $(check_exe $i) ]]; then
		msg "command '$i' can't be found"
		exit 3;
	fi
done

if [[ $# -lt 3 ]]; then
	usage;
	exit 1;
fi

keepAvg="";
noZero="";
posArgs=();
pseudoCnt=0;

while [[ $# -gt 0 ]]; 
do
	k=$1;
	shift;
	case $k in
		--keep-avg)
			keepAvg="T";
			;;
		--no-zero)
			noZero="--no-zero";
			;;
		--pseudo-count)
			pseudoCnt=$1;
			shift;
			;;
		*)
			posArgs+=("$k")
			;;
	esac
done

if [[ ${#posArgs[@]} -lt 3 ]]; then
	msg "Two reads files and one/two bed files are needed for this program"
	exit 2;
fi

set -- ${posArgs[@]}

sigFile1=$1;
sigFile2=$2;
bedFile1=$3;
bedFile2=${4:-$bedFile1};

#echo $bedFile2; exit 1;

tmpFiles=()
out=tmp.cor.$$
bed1=$out.1.bed
bed2=$out.2.bed
tmpFiles+=($bed1 $bed2)
less $bedFile1 | cut -f 1-4 >$bed1
less $bedFile2 | cut -f 1-4 >$bed2

# if bam file, convert it to bed file first
if [[ $sigFile1 =~ \.bam$ ]]; then
	tmpFile=$(rand_str).bed.gz
	tmpFiles+=($tmpFile)
	msg "Converting $sigFile1 from bam to bed [$tmpFile]"
	bam2bed.sh $sigFile1 | sort -k1,1 -k2,2n | gzip -c >$tmpFile
	sigFile1=$tmpFile
fi
if [[ $sigFile2 =~ \.bam$ ]]; then
	tmpFile=$(rand_str).bed.gz
	tmpFiles+=($tmpFile)
	msg "Converting $sigFile2 from bam to bed [$tmpFile]"
	bam2bed.sh $sigFile2 | sort -k1,1 -k2,2n | gzip -c >$tmpFile
	sigFile2=$tmpFile
fi

nf1=$(file_nf $bed1)
nf2=$(file_nf $bed2)
ok=$(add_signal $bed1 $sigFile1);
ok=$(add_signal $bed2 $sigFile2);
correlation_coef.pl $noZero --log2 $pseudoCnt --col $(( nf1 + 1)),$((nf2+1)) $bed1 $bed2

if [[ $keepAvg ]]; then
		paste $bed1 $bed2 | gzip -c \
		>$out.depth_cor.tsv.gz
		msg "Read depths for regions are saved in '$out.depth_cor.tsv.gz'"
fi

# clean up
if [[ $tmpFiles ]]; then
	rm ${tmpFiles[@]}
	#msg "No clean"
fi

msg "[`basename $0`] Job is done at `date`"

exit 0;

