#!/usr/bin/env python

import sys
import argparse as ap
import logging
import gzip

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

strands = {
        '+'  : 1,
        '+1' : 1,
        '-'  : -1,
        '-1' : -1
        }

# functions
def open_file(f, mode="r"):
    if f.endswith(".gz"):
        if 't' not in mode:
            mode += 't'
        fh = gzip.open(f,mode)
    else:
        fh = open(f,mode)
    return(fh)

desc="""
This program reads a bed file and output promoters defined by users,
i.e., the regions surrounding the start position.

Note that the promoter region will be trucated to the 3'end of a gene
if the user-specified coordinate goes beyond that.

Options (default values in []):
"""

epilog = f"""
Example:
{sys.argv[0]} --bs 1000 --as 50 -o promoters.bed.gz genes.bed.gz
"""

op = ap.ArgumentParser(
        description=desc,
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=epilog
        )

op.add_argument(
        "bedfile",
        help="input file in (gzipped) bed format"
        )

op.add_argument(
        "--bs",
        help="the number of bases before start position [%(default)d]",
        type=int,
        default=1000
        )

op.add_argument(
        "--as",
        help="the number of bases after start position [%(default)d]",
        type=int,
        dest="aS",
        default=500
        )

op.add_argument(
        "--ignore-strand",
        help="if provided, all input intervals are assumed from plus strand",
        action="store_true",
        dest="ignoreStrand"
        )

op.add_argument(
        "--skip-unstrand",
        help="if provided, intervals without strand info are skipped, otherwise plus strand is assumed",
        action="store_true",
        dest="skipUnstrand"
        )

op.add_argument(
        "-o", "--outfile",
        help="output filename [screen]",
        dest="outFile",
        default=None
        )

args = op.parse_args()

if args.outFile is None:
    outFh = sys.stdout
else:
    outFh = open_file(args.outFile,"w")

if args.ignoreStrand:
    logger.warn("Strand information for inputs is ignored")

skipped = 0
counter = 0
unstranded = 0

with open_file(args.bedfile) as fh:
    for row in fh:
        counter += 1
        if counter % 100000 == 0:
            logger.info(f"{counter} lines processed")
        fields = row.rstrip().split("\t")
        # check strands
        if args.ignoreStrand:
            strand = 1
        else:
            if len(fields) < 5:
                strand = None
            else:
                strand = fields[4]
                strand = strands[strand] if strand in strands else None
            if strand is None:
                if args.skipUnstrand:
                    skipped += 1
                    continue
                else:
                    strand = 1  # assume plus strand
                    unstranded += 1
        if strand > 0:  # plus strand
            start = int(fields[1])
            newS = start - args.bs
            newS = 0 if newS < 0 else newS
            newE = start + args.aS
            # keep newE within the gene
            newE = fields[2] if newE > int(fields[2]) else newE
        else:
            start = int(fields[2])
            newS = start - args.aS
            newS = fields[1] if newS < int(fields[1]) else newS
            newE = start + args.bs
            # Todo: need check if newE beyond chromesome length
        fields[1] = str(newS)
        fields[2] = str(newE)
        outFh.write("\t".join(fields)+"\n")

if outFh != sys.stdout:
    outFh.close()

logger.info("Job done: %d processed, %d unstranded, %d skipped" %
        (counter, unstranded, skipped))

sys.exit(0)

