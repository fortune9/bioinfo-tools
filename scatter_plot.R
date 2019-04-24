#!/bin/env Rscript

# parse configuration parameters in a file
parse_config<-function(f)
{
	# do nothing this moment
}

cur_exe<-function()
{
	args<-commandArgs(F)
	pat<-"^--file="
	path<-sub(pat,'',grep(pat,args,value=T))
}

usage<-function()
{
	exe<-cur_exe();
	message("
Usage: ", exe, " <input-file> <out-file> <x-var> <y-var> [<group> <config-file>]

This program reads an input file and make scatter plot using R.

<input-file>: the file containing the data.
<out-file>: the file to store the plots.
<x-var>: the x variable in plot, given by a column name.
<y-var>: the y variable in plot, given by a column name.
<group>: an optional argument, specifying a column which contains
the group variable. If specified, a plot for each group is also made.
Specify \"\" if no group variable.
<config-file>: a optional file with further configuration 
used for making plots. One example file will be like this:
xlim	1,10
ylim	100,200
col	blue

Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com
");
}


args<-commandArgs(T)

if(length(args) < 4)
{
	usage();
	q("no")
}

inFile<-args[1];
outFile<-args[2];
xVar<-args[3];
yVar<-args[4];
group<-args[5];
configFile<-args[6];
if(grepl("^\\s*$",group))
{
	group<-NA;
}

params<-list();
if(!is.na(configFile))
{
	params<-parse_config(configFile)
}

# global variables
inSep<-"\t";

message("Reading data")
dat<-read.delim(inFile, sep=inSep)

message("Generating plots")
ncolPlot<-2
nrowPlot<-1
if(!is.na(group))
{
	groupVals<-unique(dat[,group])
	nrowPlot<-(length(groupVals) + 1)/ncolPlot;
}

figWid<-7;
figHei<-3*nrowPlot;
#save.image("tmp.Rdata")

pdf(file=outFile, width=figWid, height=figHei);
with(dat, plot(get(xVar),get(yVar)))
dev.off()

message("Job is done")

