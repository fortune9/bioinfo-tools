#!/bin/env python

import sys;
import argparse as ap;
import textwrap;
#import primer3;
import melting;

# define functions
def calc_Tm(s, digits=3):
	''' calculate Tm using the method in melting,
	as the primer3.calcTm produces -999999.9999 for some sequences
	'''
	tm = melting.temp(s, DNA_c=50, Na_c=50, Mg_c=0, dNTPs_c=0.8);
	return(round(tm,digits));

def count_CpGs(s):
	''' count CpGs in a sequence '''
	n = s.count('CG');
	return(n);

def calc_GC_content(s):
	''' calculate GC content in a sequence '''
	n = s.count('C') + s.count('G')
	p = int(n/len(s)*100);
	return(p)

desc = textwrap.dedent("""
	This program calculates the following features for each
	sequence (unless turned off by options):

	Tm:    melting temperature
	nCpG:  number of CpGs
	GC:    the percentage of GC content
	""");

authorInfo = '''
Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

# set up arguments
optParser = ap.ArgumentParser(
		description=desc, 
		formatter_class=ap.RawTextHelpFormatter,
		#formatter_class=ap.ArgumentDefaultsHelpFormatter,
		epilog=authorInfo
		);
## positional required arguments
optParser.add_argument("infile",
	help="input file with sequences");

## optional auxillary arguments
optParser.add_argument("--column", "-c",
	help="which column contains sequences,\n default is first column, i.e., 0",
	type=int,
	metavar="column",
	dest="seqCol", # the attribute name to store value
	default=0,
	action='store',
	required = False # just for showcase purpose
	);
optParser.add_argument("--outfile", "-o",
		help="output filename",
		dest="outFile", # for demonstration only
		default=sys.stdout,
		metavar="stdout")
optParser.add_argument("--sep", "-s",
		help="field separator of input file",
		metavar="field-sep",
		default="\t"
		);
optParser.add_argument("--no-tm", "-nt",
		help="don't compute Tm",
		action="store_true") # default is false
optParser.add_argument("--no-ncpg", "-nn",
		help="don't compute nCpG",
		action="store_true") # default is false
optParser.add_argument("--no-gc", "-ng",
		help="don't compute GC content",
		action="store_true") # default is false

args = optParser.parse_args();

#print(args);

# start the analyses
o=open(args.outFile, "w")
f=open(args.infile, "r")
## determine what features to calculate
featureFuncs = {
		'Tm'   : calc_Tm,
		'nCpG' : count_CpGs,
		'GC'   : calc_GC_content


		};
features=[];
funcs=[];
for feat, func in featureFuncs.items():
	if getattr(args, 'no_' + feat.lower()):
		continue
	features.append(feat);
	funcs.append(func);

o.write(args.sep.join(["id",'seqLen'] + features)+'\n');
id=0;
for r in f:
	r = r.strip().split(args.sep);
	seq = r[args.seqCol].upper(); # uppercase sequences
#	print(seq);
	id += 1;
	res=[str(id), str(len(seq))];
	for func in funcs:
		res.append(str(func(seq)));
	#r.extend(res); # combine original line and results
	#print(args.sep.join(res))
	o.write(args.sep.join(res)+"\n");
	if id % 10000 == 0:
		print("[Info] {} sequences have been processed",
				file=sys.stderr)

f.close();
o.close();
print("Job is done", file=sys.stderr);
sys.exit(0);
