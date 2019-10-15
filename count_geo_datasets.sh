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

source $bashFuncFile && rm $bashFuncFile;

function usage()
{
	cat <<EOF
Usage: $0 [options]

This program is used to count the number of datasets published in GEO
each year of the user-specified range for given platforms.

Options (default values are in []):

-h|--help   print out this usage information.

Required:

-s|--start-year <int>  the starting year for retrieving data [2010].

-e|--end-year <int> the ending year for retrieving data [2019].

-p|--platforms <string> platform accessions for which the datasets
using the platforms are counted [].

Optional:

--species <string> if provided, only consider datasets from this
species. The parameter must follow the format used in NCBI ENTREZ
search. E.g.: "Homo sapiens" []

-o|--out-file <path> where to output results [/dev/stdout]

Example usages:

# search the datasets from 2018 to 2019 in the platforms GPL16304 and GPL13534
$0 -p GPL13534,GPL16304  -s 2018 -e 2019

# add a requirement for human samples only
$0 -p GPL13534,GPL16304 --species "Homo Sapiens" -s 2018 -e 2019

# write to output file instead of screen
$0 -p GPL13534,GPL16304 --species "Homo Sapiens" -s 2018 -e 2019 -o \
dataset_counts.tsv

EOF

}

depends=(esearch xtract)

for e in ${depends[@]}
do
	if [[ ! $(check_exe $e) ]]; then
		msg "Command '$e' does not exist"
		msg "Please install Linux package 'ncbi-entrez-direct'"
		exit 1;
	fi
done

# global variables
db=gds

# input arguments and default values
startYear=2010
endYear=2019
platforms=""
outFile=/dev/stdout
species="";


while [[ $# -gt 0 ]];
do
	k=$1; shift;
	case $k in
		-h|--help)
			usage;
			exit 3;
			;;
		-s|--start-year)
			startYear=$1;
			shift;
			;;
		-e|--end-year)
			endYear=$1;
			shift;
			;;
		-p|--platforms)
			platforms=$1;
			shift;
			;;
		--species)
			species=$1;
			shift;
			;;
		-o|--out-file)
			outFile=$1;
			;;
		*)
			otherArgs+=("$k");
			;;
	esac
done


if [[ ! "$platforms" ]]; then
	msg "The argument '--platforms' is missing"
	exit 3;
fi

# processing the arguments
platforms=($(str_split ',' $platforms))

# construct the searching strings
queryStr="GSE[etyp]"
if [[ $species ]]; then queryStr+=" AND \"$species\"[orgn]"; fi
platStr="";
for plat in ${platforms[@]}
do
	if [[ $platStr ]]; then
		platStr+=" OR \"$plat\"[ACCN]"
	else
		platStr="(\"$plat\"[ACCN]"
	fi
done

platStr+=")"; # closing parenthesis

#echo $platStr

queryStr+=" AND $platStr"

#msg "Constructed query string is: '$queryStr'"

#touch $outFile
echo -e "#year\tGSE\tsampleCount" >$outFile

for yr in `seq $startYear $endYear`;
do
	query="$queryStr AND $yr[pdat]"
	msg "Searching '$query' in database $db"
	esearch -db $db -query "$query" | esummary -format docsum \
		| xtract -pattern DocumentSummary -pfx "$yr\t" \
		-element DocumentSummary/Accession -block Samples \
		-element "#Accession" >>$outFile
		#| xtract -pattern ENTREZ_DIRECT \
		#-pfx "$yr\t" -element Count >>$outFile
	sleep 1;
done

msg "Job is done"

exit 0;


