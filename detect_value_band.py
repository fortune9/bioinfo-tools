#!/usr/bin/env python

import argparse as ap;
import sys;
import re;

# functions
def warn(s):
	''' output warning message
	'''
	print(s, file=sys.stderr);

def output(l):
	'''output a list
	'''
	print(args.sep.join(map(str, l)), file=o);

def check_len(l, size):
	'''check whether the size of provided list match 'size'
	'''
	return True if len(l) == size else False;

def clear_bucket():
	'''output the stored segments so far and reset the bucket
	'''
	global segBucket;
	if len(segBucket) < 1 or len(segBucket[0]) < 1: # no data
		return None;
	# calculate mean and stdev
	dataWidth=len(segBucket)-2; # last 2 items are lengths and range
	lengths=segBucket[len(segBucket)-2];
	means=[];
	sdevs=[];
	if args.weighted and lastChrom is not None:
		n=sum(lengths);
		for i in range(dataWidth):
			total=0;
			seg=segBucket[i];
			total=sum(seg[j]*lengths[j] for j in range(len(seg)));
			m=total/n;
			sd=(sum([(seg[j]*lengths[j]-m)**2 for j in range(len(seg))])/(n-1))**0.5 if n > 1 else 0;
			means.append("{:.5g}".format(m));
			sdevs.append("{:.5g}".format(sd));
	else:
		n=len(lengths);
		for i in range(dataWidth):
			seg=segBucket[i];
			total=sum(seg);
			m=total/n;
			sd=(sum([(x-m)**2 for x in seg])/(n-1))**0.5 if n > 1 else 0;
			means.append("{:.5g}".format(m));
			sdevs.append("{:.5g}".format(sd));
	# output the values
	chrStr = [lastChrom] if lastChrom is not None else [];
	output( chrStr + segBucket[len(segBucket)-1] + means + sdevs);
	# clear bucket
	segBucket=[];


def count_passes(datum):
	'''count the number of data points meeting the criteria
	'''
	passes=0;
	if minVals is not None and maxVals is not None:
		for i in range(len(datum)):
			if minVals[i] <= datum[i] <= maxVals[i]:
				passes+=1;
	elif minVals is not None: # only on minimal values
		for i in range(len(datum)):
			if minVals[i] <= datum[i]:
				passes+=1;
	elif maxVals is not None:
		for i in range(len(datum)):
			if datum[i] <= maxVals[i]:
				passes+=1;
	# return values
	return passes;

def grow_segs(datum, start, length):
	'''grow the bucket by appending the datum after checking all
	criteria and return regions formed.
	Return values are:
	None: datum is good
	-1: datum not passed
	'''
	if count_passes(datum) < len(datum): # bad row
		clear_bucket();
		return -1;
	# the datum is good, add it bucket
	end=start + length -1;
	if len(segBucket) < 1: # no data yet
		for i in range(len(datum)):
			segBucket.append([datum[i]]); # append the value list
		# finally add lengths and the segment range
		segBucket.append([length]); # add length list
		segBucket.append([start,start + length -1]); # standard coordinates
		return None;
	# otherwise we need check existing segs
	for i in range(len(datum)):
		segBucket[i].append(datum[i]);
	# update range values and add length
	segBucket[len(segBucket)-2].append(length);
	segRange=segBucket[len(segBucket)-1]; # last element
	if end > segRange[1]:
		segRange[1]=end; # it changes segBucket value too
	return None;


# set up arguments
desc = """
This program scans the values in the file, gets the bands of
continuous values satisfying specified criteria (see options),
and calculates some statistics for each band (currently mean
and standard deviation).

The values in multiple fields can be scanned and each field can have
its own criteria for filtering; only lines satisfying the criteria of
all specified fields are included in a band.

Default value for each option is in [];
""";

authorInfo = '''
Author:  Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

op = ap.ArgumentParser(
		description=desc,
		formatter_class=ap.RawTextHelpFormatter,
		epilog=authorInfo
		);

op.add_argument("inFile",
		help="input file"
		);

op.add_argument("--sep",
		help="field separator for input file [<tab>]",
		default="\t"
		);


op.add_argument("-f", "--fields",
		help="the field numbers from which the values will be checked. First field number is 0, second is 1, etc. [0]",
		nargs="*",
		type=int,
		dest="fields",
		default=[0],
		action='store'
		);

op.add_argument("--names",
		help="field names for the above specified fields [field1 field2 ...]",
		nargs="*",
		type=str,
		dest="fieldNames"
		)

op.add_argument("-n", "--min",
		help="the minimums (inclusive) for the examined value, in the same order as specified for --fields",
		nargs="+",
		type=float,
		dest='minCutoff'
		);

op.add_argument("-x", "--max",
		help="the maximums (inclusive) for the examined value, in the same order as specified for --fields",
		nargs='+',
		type=float,
		dest='maxCutoff'
		);

op.add_argument("-w", "--weighted",
		help="if provided, mean and deviation calculations are weighted by length (see option --start and --end). Default is false",
		action='store_true',
		dest='weighted'
		);

op.add_argument("--chr",
		help="the field number from which chromosome name is obtained. This is needed when scanning genomic regions, and input lines should be coordinate-sorted. [None]",
		default=None,
		dest='chrField',
		type=int
		);

op.add_argument("--start",
		help="the field number of genomic start position. [None]",
		default=None,
		dest='startField',
		type=int
		);
		
op.add_argument("--end",
		help="the field number of genomic end position. [None]",
		default=None,
		dest='endField',
		type=int
		);

op.add_argument("--size",
		help="the maximum size for output regions. For genomic regions, its is the length in bps; for others, it is the number of input lines. Default is no limit. [None]",
		default=None,
		dest='maxSize',
		type=float
		);

args = op.parse_args();
minVals=args.minCutoff;
maxVals=args.maxCutoff;
if minVals is None and maxVals is None:
	warn("At least one of the arguments '--min','--max' is needed");
	sys.exit(2);

if not ( all(x is None for x in	[args.chrField,args.startField,args.endField]) 
		or all(x is not None for x in [args.chrField,args.startField,args.endField])):
	warn("The arguments '--chr', '--start', and '--end' are all needed if any is provided")
	sys.exit(3);

fieldNames=args.fieldNames;

if fieldNames is None:
	fieldNames = [ 'field.' + str(i) for i in range(1, len(args.fields)+1)]

# check all the provided parameters consistent
fieldNum=len(args.fields);
for v in [minVals, maxVals, fieldNames]:
	if v is None: continue;
	if not check_len(v, fieldNum):
		warn("The size of parameter ["+' '.join(map(str, v))
				+"] does not match fields size "+str(fieldNum) );
		sys.exit(3);

i=open(args.inFile, "r");
o=sys.stdout;
header=[] if args.chrField is None else ['chrom'];
header.extend(['start','end']);
header.extend([x + '.mean' for x in fieldNames]);
header.extend([x + '.stdev' for x in fieldNames]);
output(header);

lastChrom = None;
segBucket=[]; # stores the grown segments, structure is seg1, seg2,..., segm, length-list, seg-range
lineNum=0;
lastStart=None;
for r in i:
	if re.match("^\s*$", r): # empty line
		continue;
	lineNum += 1;
	fields=r.strip().split(args.sep);
	datum=[float(fields[i]) for i in args.fields]; # get data from this line
	start=lineNum; # for non-genomic case
	datumLen=1;
	if args.chrField is not None: # genomic regions
		chrom,start,end=[fields[i] for i in [args.chrField,
			args.startField, args.endField]];
		start=int(start); # get coordinates
		end=int(end);
		if lastStart is not None and lastStart > start:
			sys.exit(
				"Input file is not coordinate sorted [{0}:{1}-{2}]".
				format(lastChrom,start,end));
		if lastChrom is not None and lastChrom != chrom: # new chrom
			clear_bucket(); #output all remaining segments
			lastChrom=chrom;
		datumLen=end-start+1;
		lastStart=start;
	# add this datum
	ok=grow_segs(datum, start, datumLen);
	# something is wrong if ok is not none
	if lineNum % 10000 == 0:
		warn(str(lineNum) + " lines have been processed");

i.close();

clear_bucket(); # don't forget to clean up

print("Job is done", file=sys.stderr);
sys.exit(0);

