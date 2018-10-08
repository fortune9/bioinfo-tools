#!/usr/bin/env python

import sys;
import os;
import argparse as ap;
import textwrap;
import pandas as pd;
import random as rd;
import bisect as bs;

# define functions
## build the map of divided regions along chromosomes
def define_regions(chrSize, trim5=0, trim3=0):
	''' return a data frame containing the starting coordinates,
	number of regions on each chromosome
	'''
	res=chrSize.apply(lambda x: [x[0],trim5,int((x[1]-trim5-trim3)/wS)], 
			axis=1,	result_type="expand");
	res.rename(columns={0: 'chr', 1:'start', 2:'counts'},
			inplace=True)
	return(res)

## get chromosome coordinates for a given region number
def region_coord(regions, number):
	'''
	give a region number, chromosome coordinates for the region are
	returned. Note the row order of the input data.frame "regions"
	matters for locating the region.
	The returned coordinates are 0-based and right open, like UCSC
	'''
	regionBounds=regions['counts'].cumsum().tolist();
	chrIndex=bs.bisect_left(regionBounds, number);
	chrInfo=regions.iloc[chrIndex,].tolist();
	preRegions=0; # the number of regions before this index
	if chrIndex > 0:
		preRegions=regionBounds[chrIndex-1];
	# now calculate chr start/end
	end=(number - preRegions)*wS + chrInfo[1];
	start=end-wS;
	return([chrInfo[0], start, end]);

def region_overlap(r1, r2):
	'''
	Get the overlap between two regions, in the format of 
	[chr, start, end, others]. The 'others' columns will be simply
	copied to results.
	'''
	if r1[0] != r2[0]:
		res=["",-1,-1,0];
	else:
		start=max(r1[1],r2[1]);
		end  =min(r1[2],r2[2]);
		overlap =end-start;
		if overlap < 0: # no overlap
			overlap=0;
			start=end=-1;
		res=[r1[0],start,end,overlap];
	res.extend(r1); # append the original records
	res.extend(r2);
	return(res);

## get the blocks overlapped with a given region
def find_overlap_blocks(blocks, region,avg=False):
	'''
	the "region" contains the coordinate, and blocks overlapped with
	the region are returned.
	'''
	subBlocks=blocks[blocks['chr']==region[0]]; # blocks on the chromosome only
	if subBlocks.empty:
		return(None);
	startIndex=bs.bisect_left(subBlocks['start'].tolist(), region[1]); # first start >= region start
	endIndex=bs.bisect_left(subBlocks['end'].tolist(), region[2]); # ends
	endIndex+=1;
	if startIndex > 0:
		startIndex-=1; # ensure all overlap blocks are included
	# get the overlapped regions
	overlappedBlocks=subBlocks.iloc[startIndex:endIndex,].apply(
			region_overlap, axis=1, result_type="expand",r2=region);
	if avg: # return the average density
		depths=overlappedBlocks.iloc[:,3]*overlappedBlocks.iloc[:,7];
		avg=depths.sum()/(region[2]-region[1]); # python 3 assumed
		return(avg)
	else: # return the blocks
		return(overlappedBlocks)
	
## get p value for a depth value
def p_for_depth(bgDepths, depth):
	'''
	Given a depth value, the fraction of values in bgDepths which are
	greater than depth is returned. bgDpeth has been sorted, so that
	we can use bisect algorithm
	'''
	if depth==None: return(None);
	index=bs.bisect_left(bgDepths,depth);
	size=len(bgDepths);
	p=(size-index)/size; # upside fraction
	if size == index: p = 1/size;
	return(p);

# set up arguments
desc = """
This program calls peaks/signal clusters from the read depth data.
It achieves the goal through two steps: first, it detects candidate
regions by sliding windows and comparing the mean read depth to the
background; second, it sharpens the ends of the regions and connects
contiguous regions.

Default optional values are in [].
""";

authorInfo = '''
Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

optParser = ap.ArgumentParser(
	description=desc,
	formatter_class=ap.RawTextHelpFormatter,
	epilog=authorInfo
		);

## positional required arguments
optParser.add_argument("infile",
		help="input file with sequence depth in bedGraph format, must be already sorted in coordinates");
optParser.add_argument("csfile",
		help="a file containing chromosome sizes in the format 'chr<tab>size' in each row. Chromosomes present in both input files will be analyzed");

## optional auxillary arguments
optParser.add_argument("-w", "--window", 
		help="the size of the window used for scanning candidate regions [300]",
		type=int,
		dest="windowSize",
		default=300,
		action='store'
		);

optParser.add_argument("-n",
		help="the number of regions sampled from genome for estimating background distribution [10000]",
		type=int,
		dest="sampleSize",
		default=10000,
		action='store'
		);

optParser.add_argument("-q","--fdr",
		help="the false discovery rate for identified peaks [0.01]",
		type=float,
		dest="qCutoff",
		default=0.01,
		action="store"
		);

optParser.add_argument("--peak-frac",
		help="the expected fraction of genome having peaks. This number is used for trimming top depth regions [0.05]",
		type=float,
		dest="peakFrac",
		default=0.05,
		action="store"
		);

optParser.add_argument("--outfile", "-o",
		help="output filename [stdout]",
		dest="outFile", # for demonstration only
		default=sys.stdout,
		metavar="stdout");

args = optParser.parse_args();
wS=args.windowSize; # the length of each region for initial scanning
sampleSize=args.sampleSize; # number of regions for background calculation
fdr=args.qCutoff; # the FDR cutoff for kept regions.
peakFrac=args.peakFrac;
if peakFrac < 0 or peakFrac > 1: 
	raise ValueError("The value for peakFrac should be in [0,1]");

o=open(args.outFile, "w");

i=1;
print("Step {0:2d}: reading read depth data and chromosome sizes".format(i));

dat=pd.read_csv(args.infile, sep="\t", header=None,
		names=["chr","start","end","depth"], skiprows=1,
			nrows=20000);
chrSizes=pd.read_csv(args.csfile,
sep="\t",header=None,names=["chr","size"]);
# only consider the chromosomes existing in the data
chrSizes=chrSizes[chrSizes['chr'].isin(dat['chr'].unique())]
if chrSizes.empty:
	sys.exit("No common chromosomes found between input files");

i+=1;
print("""
Step {0:2d}: calculate background distribution with {1} sampled
regions of length {2}""".format(i, sampleSize, wS));

# get how many non-overlapped regions exist given the window size
trimSize=10000;
regions=define_regions(chrSizes,trimSize,trimSize);
regionNum=regions['counts'].sum();
# now sample regions for background calculation
sampled=rd.sample(range(regions['counts'].sum()),int(sampleSize*(1+peakFrac)));
bgRegions=[region_coord(regions, x) for x in sampled];
bgDepths=[]; # store sorted depth values
tmp=[bs.insort(bgDepths,find_overlap_blocks(dat,r, True)) for r in bgRegions];
tmp=None;
# trim the top x% values which may be from true binding regions
bgDepths=pd.Series(bgDepths);
cutoff=bgDepths.quantile(1-peakFrac);
bgDepths=bgDepths[bgDepths<=cutoff];
## what is the distribution look like: poisson or guassian?
bgFile="bg-depth."+str(os.getpid())+".csv";
bgDepths.to_csv(bgFile, index=False);
bgDepths=bgDepths.tolist();

i+=1;
# this is the slowest step
print("Step {0:2d}: scan for significant regions with FDR={1:.2f}".format(i, fdr));
counter=0;
results=[];
for j in range(1,regionNum+1):
	r=region_coord(regions,j);
	d=find_overlap_blocks(dat,r, True);
	p=p_for_depth(bgDepths,d);
	r.extend([d,p]);
	#results.append([d,p]);
	o.write("\t".join(map(str,r))+"\n");
	counter+=1;
	if counter % 10000 == 0:
		print("{:10d} regions scanned".format(counter), file=sys.stderr);

i+=1;
print("Step {0:2d}: sharpen and merge regions".format(i));

# output the results

o.close();

print("Job is done!!", file=sys.stderr);

