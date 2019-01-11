#!/usr/bin/env python

from __future__ import print_function;
import sys;
import argparse as ap;

# functions
def warn(s):
	''' print out warning message
	'''
	print(f"[Warning] {s}", file=sys.stderr);

def output(c):
	'''output the crowdedness parameters
	'''
	print(sep.join(map(str,c)), file=o);

def calc_crowdedness(focus, before, after):
	'''calculate the following crowdness parameters
	1. the distance from the closest region, -1 if no regions within
	the considered window.
	2. number of regions within the considered window.
	3. the weighted sum of nearby regions. Now simply reciprocals of
	distances, in future, kernel-smoother-like weights can be used.
	before and after are arrays of regions before and after the focus
	region in coordinates
	'''
	dist=None;
	if len(before) > 0:
		dist = focus[1] - before[len(before)-1][2]; # start - end
	if len(after) > 0:
		tmp = after[0][1] - focus[2];
		if dist is None or tmp < dist:
			dist = tmp;
	if dist is None: dist='-1'; # no nearby regions
	nr = len(before) + len(after);
	ws = sum([1.0/d if d > 0 else 1 for d in [ focus[1] - r[2] for r in before] ] +
			[1.0/d if d > 0 else 1 for d in [ r[1] - focus[2] for r in after] ]);
	return( focus + [dist, nr, ws]); # append the data

def shift_focus(focus, before, after):
	'''
	shift focus to next region, and clean up the 'before stack'
	'''
	before.append(focus); # add current focus to before regions
	if len(after) < 1: return (None,[],[]); # no more focus
	focus=after[0];
	after=after[1:];
	before=list(filter(lambda x:focus[1]-x[2]<beforeDist,
			before)); # remove out-of-range before-regions
	return (focus, before, after);

def finish_left_regions(focus, before, after):
	'''calculate crowdness for current focus and all regions in
	the 'after' stack
	'''
	res=[];
	tmp=calc_crowdedness(focus, before, after);
	output(tmp);
	res.append(tmp);
	while len(after) > 0:
		focus, before, after = shift_focus(focus, before, after);
		tmp=calc_crowdedness(focus, before, after);
		output(tmp);
		res.append(tmp);
	return(res);

desc='''
Given a list of genomic regions, this program calculates the
parameters to measure how close each region is to nearby ones.
The parameters are:

1. minD, the distance between a considered region and the closest one.
2. nr, the number of region within a local window (see option
--window-len).
3. sr, the sum of reprocicals of distances between the considered
region and the others within a local window.

Default values for optional arguments are in [].
'''

authorInfo='''
Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

op=ap.ArgumentParser(
		description=desc,
		formatter_class=ap.RawTextHelpFormatter,
		epilog=authorInfo);

op.add_argument("inFile",
		help="the input file containing region coordinates in bed format and sorted according to coordinates");

op.add_argument("--out-file", "-o",
		help="output filename [stdout]",
		default=sys.stdout,
		dest="outFile"
		);

op.add_argument("--window-len", "-w",
		help="the size of a local window used to calculate region density; a larger size anticipates smoother values across regions. If the size is 300, then the regions within 150bp from both sides are included [300]",
		type=int,
		dest="winLen",
		default=300,
		action='store',
		required=False
		);

op.add_argument("--sep", "-s",
		help="field separate for input file [<tab>]",
		default="\t",
		action="store"
		);

args=op.parse_args();
winLen=args.winLen;
beforeDist=afterDist=winLen/2.0;
sep=args.sep;

i=open(args.inFile, "r");
o=open(args.outFile, "w");
output(['chr','start','end','name','minD','nr','sr']);

lastStart=None;
lastChr=None;
analyzedChroms=[];
focus=None;
beforeRegions=[];
afterRegions=[];
counter=0;

for r in i:
	r = r.rstrip().split(sep);
	if len(r) < 4:
		warn("Less than 4 fields found on " + sep.join(r));
		sys.exit(1);
	chrom, start, end, name = r[:4];
	if not start.isdigit():
		warn(f"This '{start}' is not number for start");
		continue;
	counter += 1;
	if counter % 2000 == 0:
		warn(f"{counter} lines have been processed");
	if lastChr is not None and lastChr != chrom: # new chromosome
		if chrom in analyzedChroms:
			sys.exit(f"chromosome '{chrom}' appears in multiple blocks. file not sorted")
		analyzedChroms.append(chrom);
		# This is new chromosome and we need finish all the remaining
		# regions
		resSet=finish_left_regions(focus, beforeRegions, afterRegions);
		# reset all regions
		focus=None;
		beforeRegions=[];
		afterRegions=[];
	
	start = int(start);
	end = int(end);
	if focus is None: # no focused region yet
		focus=[chrom,start,end,name];
		lastChr=chrom;
		lastStart=start;
		continue;
	# there is focus, check the distance to it
	if lastChr == chrom and start < lastStart:
		sys.exit("The input file is not coordinate sorted");
	## the new region is out of range, run calculation for current
	## focused region and shift to next region until the new region
	## is within distance
	while start-focus[2] >= afterDist:
		res=calc_crowdedness(focus, beforeRegions, afterRegions);
		output(res);
		focus, beforeRegions, afterRegions = shift_focus(focus,
				beforeRegions, afterRegions);
		if focus is None: # no more regions in after
			break
	
	if focus is None:
		focus = [chrom,start,end,name];
	else: ## otherwise append this region to afterRegions stack only
		afterRegions.append([chrom,start,end,name]);
	lastChr=chrom;
	lastStart=start;

# finish the remaining regions
if focus is not None: finish_left_regions(focus, beforeRegions, afterRegions);

i.close();
o.close();

warn("Job is done");

sys.exit(0);

