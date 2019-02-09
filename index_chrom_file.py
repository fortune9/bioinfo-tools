#!/usr/bin/env python

import argparse as ap;
import sys;

# set up arguments
desc='''
This program accepts a bed file and generates a text file showing the
byte offset for each chromosome in the bed file.

default values for arguments are shown in [].
''';

author='''
Author: Zhenguo Zhang
Email:  zhangz.sci@gmail.com
''';

op=ap.ArgumentParser(
		description=desc,
		epilog=author,
		formatter_class=ap.RawTextHelpFormatter
		);

op.add_argument("bedFile",
		help="the input file to be indexed"
		);

op.add_argument("-o", "--outfile",
		help="the file to store output [stdout]",
		default=sys.stdout,
		dest='outFile'
		);

args=op.parse_args();
o=args.outFile;
if type(o) is str:
	o=open(o, "w");

lastChrom=None;
chrStart=-1;
chrEnd=-1;
chrOffset=-1;
print("#" + " ".join(sys.argv), file=o);
print("\t".join(["chrom","start","end","offset"]), file=o);
lineNum=0;
with open(args.bedFile, 'r') as f:
	while True:
		r=f.readline();
		if r == '': break; # end of file
		chrom,start,end=r.rstrip().split("\t")[:3];
		if chrom != lastChrom: # new chromosome
			if lastChrom is not None: # old chromosome's record
				print("\t".join(map(str, 
					[lastChrom, chrStart, chrEnd,chrOffset])), file=o);
			# initialize for new chromosome
			chrOffset=f.tell()-len(r.encode('utf-8'));
			lastChrom=chrom;
			chrStart=start;
			chrEnd=end;
		else: # the same chrom
			chrEnd=end; # update end 
		lineNum+=1;
		if lineNum % 100000 == 0:
			print(f"[INFO] {lineNum} lines have been processed",
					file=sys.stderr);

# don't forget the last chrom
if lastChrom is not None:
	print("\t".join(map(str,[lastChrom, chrStart, chrEnd,chrOffset])),
		file=o);

if o is not sys.stdout: o.close();

print("Job done", file=sys.stderr);

sys.exit(0);

