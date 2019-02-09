#!/usr/bin/env python

import sys;
import argparse as ap;
import subprocess as sp;

# functions
def CpG_pos(s):
	'''
	find all CpGs in a sequence
	'''
	last_found=-1;
	s=s.upper();
	while True:
		last_found=s.find('CG', last_found+1);
		if last_found == -1:
			break; # no more occurences
		yield last_found;


def get_seq(chrom, start, end):
	'''
	Get sequence for a given region
	'''
	seq=sp.run(['twoBitToFa', args.twobitFile, 'stdout',
		"-seq="+chrom,"-start="+start,"-end="+end], stdout=sp.PIPE,
		universal_newlines=True);
	if seq.stdout == '': return(None);
	seq=seq.stdout.split("\n")[1:];
	seq=''.join(seq);
	return(seq);

desc='''
This program finds all the CpG positions in given genomic regions.
The output is in bed format.
''';

authorInfo='''
Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

# set up arguments
optParser=ap.ArgumentParser(
		description=desc,
		formatter_class=ap.RawTextHelpFormatter,
		epilog=authorInfo);

## positional arguments
optParser.add_argument("regionFile",
		help="file containing genomic region in bed format, having at least 4 fields");

optParser.add_argument("twobitFile",
		help="the .2bit file from which sequences can be extracted");

optParser.add_argument("-c", "--count-only",
		help="if provided, only the number of CpGs in each region is returned",
		action='store_true',
		dest='countOnly');

args=optParser.parse_args();
sep="\t";

r=open(args.regionFile,"r");
counter=0;
for l in r:
	counter+=1;
	if counter % 1000 == 0:
		print("# Start processing line {0}".format(counter),
				file=sys.stderr);
	l=l.strip().split(sep);
	(chrom, start, end, name)=l[:4];
	seq=get_seq(chrom, start, end);
	if seq == None:
		print("Cannot get sequence for [{0}]".format(l), file=sys.stderr);
		next;
	poses=list(CpG_pos(seq));
	if args.countOnly: # count only 
		print(sep.join(map(str, [chrom, start,end,name,len(poses)])))
		continue;
	# report positions otherwise
	if len(poses) < 1: # no CpGs
		print('#'+sep.join(map(str, [chrom, start,
				end,name+".CpG.NA"])));
		continue;
	else:
		start=int(start);
		i=0;
		for p in poses:
			i+=1;
			print(sep.join(map(str, [chrom, start+p,
				start+p+2,name+".CpG."+str(i)])));
r.close()

print("Job is done", file=sys.stderr);

sys.exit(0);

