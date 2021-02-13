#!/usr/bin/env python

import sys
import argparse as ap
import logging
import gzip

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


# functions
def read_filter_vals(f):
    fh = open_file(f)
    d = dict()
    for l in fh:
        d[l.rstrip()] = 1;

    return(d)

def open_file(f, mode="rt"):
    if f.endswith(".gz"):
        fh = gzip.open(f,mode)
    else:
        fh = open(f,mode)
    return(fh)

desc = """
This program filters input file based on the provided filtering values
"""

op = ap.ArgumentParser(
        description=desc,
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog=""
        )

op.add_argument(
        "infile",
        help="input file to which filter will be applied"
        )

op.add_argument(
        "-o", "--outfile",
        help="output file name. Default is screen",
        dest="outFile"
        )

op.add_argument(
        "--val-file",
        help="file containing the values, one value per line",
        dest="valFile"
        )

op.add_argument(
        "-c, --column",
        help="""the column number for the input file where the compared
        values are located [%(default)d]""",
        type=int,
        dest="valCol",
        default=1
        )

args = op.parse_args()

if args.valFile is None:
    logger.error("The option --val-file is required")
    raise SystemExit

if args.outFile is None:
    writer = sys.stdout
else:
    writer = open_file(args.outFile, "wt")

valDict = read_filter_vals(args.valFile)
valCol = args.valCol - 1

fh = open_file(args.infile)

for row in fh:
    fields = row.rstrip().split("\t")
    if fields[valCol] in valDict:
        writer.write(row)

fh.close()

logger.info("Job is done")

sys.exit(0)

