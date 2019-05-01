#!/usr/bin/env python

import sys;
import argparse as ap;


desc='''
This program calculates the group weighted mean by going through the
lines of the input file. Use the options '--group', '--weight', and
'--value' to specify which columns contain the group-id, weights, and
values, respectively.

Default optional values are in [].
''';

authorInfo = '''
Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com
''';

optParser=ap.ArgumentParser(
		description=desc,
		formatter_class=ap.RawTextHelpFormatter,
		epilog=authorInfo);

optParser.add_argument("infile",
		help="input file, the same group must be in contiguous lines"
		);

## mandatory options
optParser.add_argument("--group","-g",
		help="the column containg group id, the lines within the same group are used for calculating average, the first column is at 1",
		type=int,
		metavar="group-id",
		action="store",
		required=True);

optParser.add_argument("--weight","-w",
		help="the column containg weights",
		type=int,
		metavar="weight",
		action="store",
		required=True);

optParser.add_argument("--value","-v",
		help="the column containg values to compute average",
		type=int,
		metavar="value",
		action="store",
		required=True);

## auxilary options
optParser.add_argument("--sep", "-s",
		help="field separator of input file [tab]",
		metavar="field-sep",
		default="\t");

optParser.add_argument("--outfile", "-o",
		help="output filename [stdout]",
		dest="outFile", # for demonstration only
		default=sys.stdout,
		metavar="output");

args=optParser.parse_args();
gCol=args.group - 1;
vCol=args.value -1;
wCol=args.weight - 1;

# start analysis
o=open(args.outFile, "w") if args.outFile != sys.stdout else args.outFile;
i=open(args.infile,"r");
## read each group and calculates the average
lastG="";
totalW=0;
totalV=0;
counter=0;
for r in i:
	r=r.strip().split(args.sep);
	g=r[gCol];
	v=float(r[vCol]);
	w=float(r[wCol]);
	if lastG != g: # a new group
		if lastG != "":
			out=args.sep.join([lastG, "{:.4g}".format(totalV/(totalW or 1))]);
			o.write(out+"\n");
			counter+=1;
			if counter % 10000 ==0:
				print("{0:6d} groups have been processed\n".format(counter),
					file=sys.stderr
					);
		# initialize
		lastG=g;
		totalW=totalV=0;
	# otherwise accumulate values
	totalW += w;
	totalV += w*v;

if lastG != "":
	out=args.sep.join([lastG, "{:.4g}".format(totalV/(totalW or 1))]);
	o.write(out+"\n");
	counter+=1;

o.close();
i.close();

print("Job is done [{0} groups in total]\n".format(counter),
		file=sys.stderr);
sys.exit(0);

