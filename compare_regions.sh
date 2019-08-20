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
Usage: $0 [options] <bed-file1> <bed-file2>

**Note**: the input bed files need be coordinate sorted using the
command 'sort -k1,1 -k2,2n'.

This program accepts two files in bed format, and do
the following:

1. get the number of overlapped and specific regions.
2. get the overlapped and specific intervals. This is
different from 1 as it splits a region into intervals,
which are either overlapped with the regions from the
other file, or specific to this file.
3. signal correlations if the signal is provided (see
the options below); usually signal is read depth or
alike.

Options (default values are in []):

--fracs <num1,num2>:  two numbers in the range from 
0 to 1 to indicate how much fractions in overlap for
a region to be considered overlapped with the other
region. The two numbers are for regions in <bed-file1>
and <bed-file2>, respectively. And either of the condition
being true will yield an overlap classification. [0.5,0.5]

--signals <path1[,path2]>: one or two signal files,
which can be in the format bam/bed/gff/vcf; based on
the covered regions by these files, then signal values
are computed for input regions. If two files are provided,
then the regions from <bed-file1> will use the former
file while those from <bed-file2> will use the latter.
Otherwise both use the same file. These signal files also
need be sorted based on chromosome name and then start positions
for faster computing, except for bam files, which will be sorted
by the program internally. Otherwise the program may fail []

--keep-avg: a witch option, when provided, the files containing
computed average signals for each common peak and common interval
will be kept.

--out <string>:  the prefix for output files. Default uses
current PID [bed1_vs_bed2.\$\$].

Example uses:


EOF

rm $bashFuncFile

}

depends=(regions_overlap.sh merge_regions.sh python3 
weighted_avg.py correlation_coef.pl bam2bed.sh)

for i in ${depends[@]}
do
	if [[ ! $(check_exe $i) ]]; then
		msg "command '$i' can't be found"
		exit 3;
	fi
done

if [[ $# -lt 2 ]]; then
	usage;
	exit 1;
fi

keepAvg="";
fracs="0.5,0.5"
sigFiles=""
out="bed1_vs_bed2.$$";
posArgs=();

while [[ $# -gt 0 ]]; 
do
	k=$1;
	shift;
	case $k in
		--fracs)
			fracs=$1;
			shift;
			;;
		--signals)
			sigFiles=$1;
			shift;
			;;
		--keep-avg)
			keepAvg="T";
			;;
		--out)
			out=$1;
			shift;
			;;
		*)
			posArgs+=("$k")
			;;
	esac
done

if [[ ${#posArgs[@]} -lt 2 ]]; then
	msg "Two bed files are needed for this program"
	exit 2;
fi

set -- ${posArgs[@]}

bedFile1=$1;
bedFile2=$2;
fracs=(`echo $fracs | sed -e 's/,/\n/g'`)
frac1=${fracs[0]}
frac2=${fracs[1]}

if [[ $(pass_bc "$frac1 < 0") -gt 1 ]] || 
	[[ $(pass_bc "$frac1 > 1") -gt 1 ]]; then
	msg "The first fraction '$frac1' out of range [0,1]"
	exit 3;
fi

if [[ $(pass_bc "$frac2 < 0") -gt 1 ]] || 
	[[ $(pass_bc "$frac2 > 1") -gt 1 ]]; then
	msg "The second fraction '$frac2' out of range [0,1]"
	exit 3;
fi

msg "Step 1: compute overlapped regions"
regions_overlap.sh $bedFile1 $bedFile2 $frac1 $frac2 T $out

msg "Step 2: compute overlapped intervals"
itvFile=$out.merged_intervals.bed.gz
merge_regions.sh $bedFile1 $bedFile2 | gzip -c \
>$itvFile

if [[ $sigFiles ]]; then
	msg "Step 3: compute signal correlations"
	corrFile=$out.sig_cor.txt
	echo "# Input signal files: $sigFile1 and $sigFile2" >$corrFile
	echo "# Correlation for common regions:" >>$corrFile
	sigFiles=($(str_split "," $sigFiles))
	#echo "${sigFiles[@]}"; exit 10;
	singleFile="";
	if [[ ${#sigFiles[@]} -lt 2 ]]; then
		singleFile="T"
		sigFile1=${sigFiles[0]}
		sigFile2=$sigFile1
	else
		sigFile1=${sigFiles[0]}
		sigFile2=${sigFiles[1]}
	fi
	#echo "$sigFile1 : $sigFile2"; exit 10;
	# if bam file, convert it to bed file first
	if [[ $sigFile1 =~ \.bam$ ]]; then
		tmpFile=$(rand_str).bed.gz
		tmpFiles+=($tmpFile)
		msg "Converting $sigFile1 from bam to bed [$tmpFile]"
		bam2bed.sh $sigFile1 | sort -k1,1 -k2,2n | gzip -c >$tmpFile
		sigFile1=$tmpFile
	fi
	if [[ ! $singleFile ]]; then # two signal files
		if [[ $sigFile2 =~ \.bam$ ]]; then
			tmpFile=$(rand_str).bed.gz
			tmpFiles+=($tmpFile)
			msg "Converting $sigFile2 from bam to bed [$tmpFile]"
			bam2bed.sh $sigFile2 | sort -k1,1 -k2,2n | gzip -c >>$tmpFile
			sigFile2=$tmpFile
		fi
	else
		sigFile2=$sigFile1;
	fi
	msg "Step 3.1: correlation for common regions"
	# get common regions for each input
	nf1=$(file_nf $bedFile1)
	nf2=$(file_nf $bedFile2)
	commBedFile=$out.common.bed
	commBedFile1=$out.common.$$.1.bed
	commBedFile2=$out.common.$$.2.bed
	ef1=4; # the last field to extract
	if [[ $ef1 -gt $nf1 ]]; then
		ef1=$nf1;
	fi
	cut -f 1-$ef1 $commBedFile >$commBedFile1
	sf=$(( nf1 + 1 ));
	ef2=4;
	if [[ $ef2 -gt $nf2 ]]; then
		ef2=$nf2;
	fi
	ef2=$(( ef2 + nf1 ))
	cut -f $sf-$ef2 $commBedFile >$commBedFile2
	ok=$(add_signal $commBedFile1 $sigFile1);
	ok=$(add_signal $commBedFile2 $sigFile2);
	correlation_coef.pl --col $(( ef1 + 1)),$((ef2-nf1+1)) \
	$commBedFile1 $commBedFile2 >>$corrFile
	if [[ $keepAvg ]]; then
		paste $commBedFile1 $commBedFile2 | gzip -c \
		>$out.common.avg_sig.tsv.gz
		msg "Average signals for common peaks are in '$out.common.avg_sig.tsv.gz'"
	fi
	rm $commBedFile1 $commBedFile2
	msg "Step 3.2: correlation for common intervals"
	if [[ $singleFile ]]; then
		msg "Skipping step 3.2 as only one signal file is provided"
	else
		# get common intervals
		commItvFile=$out.common_itv.avg_sig.tsv
		less $itvFile | gawk -v FS="\t" '$4 ~ /^comm\./' >$commItvFile
		nf=$(file_nf $commItvFile)
		ok=$(add_signal $commItvFile $sigFile1);
		ok=$(add_signal $commItvFile $sigFile2);
		echo "# correlation for common intervals" >>$corrFile
		correlation_coef.pl --col $((nf+1)),$((nf+2)) \
			$commItvFile >>$corrFile
		if [[ $keepAvg ]]; then
			msg "Average signals for common intervals are in '$commItvFile.gz'"
			gzip -S .gz $commItvFile
		else
			tmpFiles+=($commItvFile)
		fi
	fi
	msg "Signal correlation coefficents are written into $corrFile"
fi

# clean up
if [[ $tmpFiles ]]; then
	rm ${tmpFiles[@]}
	#msg "No clean"
fi

msg "Job is done at `date`"

exit 0;

