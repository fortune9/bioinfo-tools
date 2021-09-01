#!/bin/bash

set -e

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <region-file> <out-file> [<pos-only>]"
    exit 1;
fi

inFile=$1
outFile=$2
posOnly=$3

posCol=""
if [[ "$posOnly" ]]; then
    posCol="--pos-col 2"
fi

echo "Adding annotation: gene"
./annotate_regions.R --out-file tmp1.$$.tsv.gz --feat-name gene \
$posCol $inFile genes.bed.gz
echo "Adding annotation: exon"
./annotate_regions.R --out-file tmp2.$$.tsv.gz --feat-name exon \
$posCol tmp1.$$.tsv.gz exons.bed.gz
echo "Adding annotation: intron"
./annotate_regions.R --out-file $outFile --feat-name intron \
    $posCol tmp2.$$.tsv.gz introns.bed.gz && \
    rm tmp1.$$.tsv.gz tmp2.$$.tsv.gz

echo "Job done"

exit 0;


