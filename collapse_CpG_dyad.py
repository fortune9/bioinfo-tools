#!/usr/bin/env python

import argparse as ap;
import sys;
import re;
import datetime as dt;

# functions
def warn(s):
	'''display warning message
	'''
	print("[WARNING]"+s, file=sys.stderr);

def output_stack(stack):
	'''collapse the rows in the stack and output the result
	'''
	if len(stack) < 1: return None; # no data
	if len(stack) == 1: # only one row
		row=stack[0];
		res=[row[i] for i in range(len(row)) 
				if i not in	[chrCol,posCol]];
	elif len(stack) == 2: # two Cs
		row1,row2 = stack;
		res=[ args.joinSep.join([row1[i],row2[i]]) for i in
			range(len(row1)) if i not in [chrCol,posCol] ];
	else:
		sys.exit(f"two many rows for position {lastChrom}:{lastPos}");
	print("\t".join(map(str, [lastChrom,lastPlusCPos] + res)),
		file=o);

def read_CpG_pos(chrom):
	''' read all CpG positions for a chromosome.
	Note: the positions are 0-based
	'''
	global CpG_poses_in_chrom;
	if not (chrom in chromOffsets):
		print(f"Chromosome '{chrom}' can't be found in CpG index file", 
			file=sys.stderr);
		sys.exit(1);
	if args.verbose:
		warn(f"Reading CpG positions for {chrom}");
	offset=chromOffsets[chrom];
	cpgPosFh.seek(offset);
	CpG_poses_in_chrom = {}; # use a hash to increase speed
	# check the first line, sometimes it can seek to wrong places
	r=cpgPosFh.readline().rstrip().split("\t");
	if len(r) < 2 or r[0] != chrom or (not re.match("^\d+$",r[1])):
		sys.exit(f"The offset '{offset} for {chrom}' seeks an unknown place " 
				+ "in cpgPosFile: " + "\t".join(r));
	CpG_poses_in_chrom[int(r[1])]=1;
	while True:
		r=cpgPosFh.readline();
		if r == '': break; # end of file
		r = r.rstrip().split("\t");
		if r[0] != chrom: break; # reach a different chrom
		CpG_poses_in_chrom[int(r[1])]=1;
	if len(CpG_poses_in_chrom)==0:
		warn(f"No CpG positions info for chromosome {chrom}");
	if args.verbose:
		warn(f"Reading CpG positions for {chrom} done");
	return 1;

def is_plus_c(chrom, pos):
	''' check whether the given position is the
	plus C in a CpG dyad
	'''
	if chrom != lastChrom: # a new chromosome
		read_CpG_pos(chrom);
	if pos in CpG_poses_in_chrom:
		return True;
	if pos-1 in CpG_poses_in_chrom:
		return False; # minus C
	return None; # not in a CpG

# set up arguments
desc="""
This program collapses two rows in input file if they represent two Cs
from the same CpG dyad.

It will report one position for each CpG dyad and join values from
different rows for each other field (see the option --join-sep for
joining separator).

Default values for optional arguments are in [].
""";

author="""
Author:  Zhenguo Zhang
Email: zhangz.sci@gmail.com
""";

op=ap.ArgumentParser(
		description=desc,
		epilog=author,
		formatter_class=ap.RawTextHelpFormatter
		);

op.add_argument("dataFile",
		help="input file containing data to collapse, must be sorted by coordinates"
		);

op.add_argument("cpgPosFile",
		help="the file containing CpG positions in a genome, in bed format"
		);

op.add_argument("chromIndexFile",
		help="file containing byte offset for each chromosome in cpgPosFile, with 1st and 2nd columns being chromosome and offsset, respectively, field separated by tab"
		);

op.add_argument("-c", "--chr",
		help="the column number for the chromosome field in dataFile. Default is 1st column [1]",
		type=int,
		default=1,
		dest="chrCol"
		);

op.add_argument("-p", "--pos",
		help="the column number for the C's genomic position in dataFile. [2]",
		type=int,
		default=2,
		dest="posCol"
		);

op.add_argument("-1","--one-based",
		help="logical option indicating whether positions in dataFile are in normal 1-based format, default is 0-based UCSC format [False]",
		action="store_true",
		dest="oneBased"
		);

op.add_argument("-s", "--sep",
		help="the field separator for dataFile. [tab]",
		default="\t",
		dest="inSep"
		);

op.add_argument("--join-sep",
		help="the separator to join values from different rows for a CpG dyad [,]",
		default=",",
		dest="joinSep"
		);

op.add_argument("-o","--outfile",
		help="the file to store output [stdout]",
		default=sys.stdout,
		dest="outFile"
		);

op.add_argument("-v", "--verbose",
		help="print more details of program running stages [False]",
		action='store_true',
		dest='verbose'
		);

args=op.parse_args();
chrCol=args.chrCol-1; # convert to python array index
posCol=args.posCol-1;

# global variables
CpG_poses_in_chrom={}; # use dict to increase search speed, faster than list
chromOffsets=dict();
# file handles
o=open(args.outFile, "r") if type(args.outFile) is str else args.outFile;
dataFh=open(args.dataFile, "r");
cpgPosFh=open(args.cpgPosFile,"r");
cpgIndexFh=open(args.chromIndexFile,"r");
for r in cpgIndexFh:
	r=r.rstrip().split("\t");
	chromOffsets[r[0]] = int(r[1]);
lineNum=0;
lastChrom="";
lastPos=-100;
lastPlusCPos=-1; # the position of plus C in last CpG
lastDyadDone=True; # indicator whether last CpG dyad completed
rowStack=[];

for r in dataFh:
	if  re.match("^\s*$", r): # empty line
		continue;
	lineNum += 1;
	fields=r.rstrip().split(args.inSep);
	chrom=fields[chrCol];
	pos=int(fields[posCol]);
	if args.oneBased: pos -= 1; # convert it to 0-based coordinates
	isPlusC=is_plus_c(chrom, pos); # test whether C is from plus strand
	if isPlusC is None: sys.exit(f"'{chrom}:{pos}' is not in CpG");
	if not isPlusC and pos - lastPos == 1 and chrom == lastChrom: #	minus C from the dyad
		rowStack.append(fields); # add it to the stack
	else: # for other new minus C, plus C
		if len(rowStack) > 0: output_stack(rowStack); # clean stack
		rowStack=[fields]; # initialize stack
		lastPlusCPos = pos if isPlusC else pos -1;
	# update position info
	lastChrom=chrom;
	lastPos=pos;
	if lineNum % 100000 == 0:
		warn(f"{lineNum} lines have been processed");

if len(rowStack) > 0: output_stack(rowStack); # last dyad

dataFh.close();
cpgIndexFh.close();
cpgPosFh.close();
if o != sys.stdout: o.close();

warn("Job is done at " + str(dt.datetime.now()));

sys.exit(0);


