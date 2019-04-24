#!/bin/env Rscript

cur_exe<-function()
{
	args<-commandArgs(F);
	path<-sub("^--file=", '', grep('^--file=', args, value=T))
}

usage<-function()
{
	exe<-cur_exe();
	message("
Usage: ", exe, " <file1> <file2> <commIds1> [<commIds2>]

This program reads two files and merge them on the common columns
specified by the arguments <commIds1> and <commIds1> (if given). 
Multiple columns can be specified and separated by comma. <commIds1>
and <commIds2> are column Ids from <file1> and <file2>, respectively.
If <commIds2> is omitted, <commIds1> is used for file <file2>.

Input <file1> and <file2> use <tab> as field separator, must have
header line, and can't be compressed.

The output is written to screen.

E.g.: ",exe," f1.bed f2.bed id11,id12 id21,id22

Author: Zhenguo Zhang
Email: zhangz.sci@gmail.com

");

}

library(data.table)

args<-commandArgs(T);
if(length(args) < 3)
{
	usage();
	q("no")
}

f1<-args[1];
f2<-args[2];
comCols1<-args[3];
comCols2<-args[4];
if(is.na(comCols2))
{
	comCols2<-comCols1;
}

comCols1<-unlist(strsplit(comCols1, ","))
comCols2<-unlist(strsplit(comCols2, ","))
if(length(comCols1) != length(comCols2))
{
	stop("The numbers of provided columns for two files differ")
}

dat1<-fread(f1)
dat2<-fread(f2)

dat<-merge(dat1, dat2, by.x=comCols1, by.y=comCols2, all=T)

fwrite(dat,file="",sep="\t")

message("Job is done")

