#!/bin/bash

set -e

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <region-file> <out-file> [<pos-only>]"
    echo "Add annotation info: gene, promoter, exon, intron, cgi"
    echo "<pos-only>: if provided, meaning <region-file> contains"
    echo "single-column position, rather than two-column bed file"
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
./annotate_regions.R --out-file tmpOut.$$.tsv.gz --feat-name gene \
$posCol $inFile genes.bed.gz && mv tmpOut.$$.tsv.gz tmpIn.$$.tsv.gz
echo "Adding annotation: promoter"
./annotate_regions.R --out-file tmpOut.$$.tsv.gz --feat-name promoter \
$posCol tmpIn.$$.tsv.gz promoters.bed.gz && mv tmpOut.$$.tsv.gz tmpIn.$$.tsv.gz
echo "Adding annotation: exon"
./annotate_regions.R --out-file tmpOut.$$.tsv.gz --feat-name exon \
$posCol tmpIn.$$.tsv.gz exons.bed.gz && mv tmpOut.$$.tsv.gz tmpIn.$$.tsv.gz
echo "Adding annotation: intron"
./annotate_regions.R --out-file tmpOut.$$.tsv.gz --feat-name intron \
$posCol tmpIn.$$.tsv.gz introns.bed.gz && mv tmpOut.$$.tsv.gz tmpIn.$$.tsv.gz
echo "Adding annotation: CGI"
./annotate_regions.R --out-file $outFile --feat-name CGI \
    $posCol tmpIn.$$.tsv.gz CGIs.bed.gz && \
    rm tmpIn.$$.tsv.gz

echo "Job done"

exit 0;


