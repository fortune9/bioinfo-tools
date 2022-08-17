#!/bin/bash

id=$1
fq11=$2
fq12=$3
fq21=$4
fq22=$5

#outDir='s3://zymo-epiquest2/zr4730/cfblood/fastq/'

if [[ $# -lt 5 ]]; then
    echo "Need 5 input parameters"
    exit 2;
fi

echo "Processing $id at `date`"

for f in $fq11 $fq21 $fq12 $fq22
do
    if [[ ! -f "$f" ]]; then
        echo "File '$f' doesn't exist"
        exit 1;
    fi
done

zcat $fq11 $fq21 | gzip -c >$id.R1.fastq.gz
zcat $fq12 $fq22 | gzip -c >$id.R2.fastq.gz

echo "Processing $id is done at `date`"

exit 0;

