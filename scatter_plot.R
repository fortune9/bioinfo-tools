#!/bin/env Rscript

# parse configuration parameters in a file
parse_config<-function(f)
{
	params<-read.table(f, sep="\t", stringsAsFactors=F)
	res<-list()
	tmp<-sapply(seq_len(nrow(params)), 
				function(i) res[params[i,1]]<<-params[i,2])
	return(res)
}

## convert a string to vector
str_to_vec<-function(s,sep=",",type=as.numeric)
{
	v<-unlist(strsplit(s,sep));
	v<-sapply(v, type);
	return(v);
}

## the function to make plot
scatter_plot<-function(x,y,statPos="topleft",...)
{
	cor.test(x, y)->corT
	r<-sprintf("%.2g", corT$est)
	p<-sprintf("%.2g", corT$p.val)
	plot(x,y,...)
	legend(statPos,legend=c(as.expression(bquote(italic(R)==.(r))),as.expression(bquote(italic(P)==.(p)))),bty="n",pch=NA,cex=0.8)
	return(NULL)
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

# the acceptable parameters in config file
knownParams<-c("xlim","ylim","xlab","ylab","main","col","figwid","plothei","ncolplot")
params<-list();
if(!is.na(configFile))
{
	params<-parse_config(configFile); # return a list
}

for(p in knownParams)
{
	assigned<-p %in% names(params);
	params[p]<-switch(p,
		   xlim=ifelse(assigned,
					   list(str_to_vec(params[[p]])),list(NULL)),
		   ylim=ifelse(assigned,
					   list(str_to_vec(params[[p]])),list(NULL)),
		   xlab=ifelse(assigned, params[p],list(NULL)),
		   ylab=ifelse(assigned, params[p],list(NULL)),
		   main=ifelse(assigned, params[p],list(NULL)),
		   col=ifelse(assigned,
					  list(str_to_vec(params[[p]],type=as.character)),"black"),
		   figwid=ifelse(assigned, as.numeric(params[[p]]),7),
		   plothei=ifelse(assigned, as.numeric(params[[p]]),3),
		   ncolplot=ifelse(assigned, as.integer(params[[p]]),2) # height of each plot
		   );
}

#save.image("tmp.Rda")
#q("no")
# global variables
inSep<-"\t";

message("Reading data")
dat<-read.delim(inFile, sep=inSep, stringsAsFactors=F)

message("Generating plots")
ncolPlot<-params$ncolplot
#message("ncolplot->", ncolPlot);
nrowPlot<-1
groupVals<-c()
if(!is.na(group))
{
	message("Group variable ", group, " is used\n")
	stopifnot(group %in% colnames(dat))
	groupVals<-unique(dat[,group])
	nrowPlot<-ceiling((length(groupVals) + 1)/ncolPlot);
}

figWid<-params$figwid;
figHei<-params$plothei*nrowPlot;
#save.image("tmp.Rdata")

pdf(file=outFile, width=figWid, height=figHei);
par(mar=c(3,3,1.1,1.1))
par(mgp=c(2, 0.8, 0))
layout(matrix(seq_len(ncolPlot*nrowPlot),nc=ncolPlot,byrow=T))
plotCols<-rep(params$col, length=ncolPlot*nrowPlot);
with(dat, scatter_plot(get(xVar),get(yVar),
	xlim=params$xlim, ylim=params$ylim,
	xlab=params$xlab, ylab=params$ylab,
	main=params$main, col=plotCols[1],
	cex=0.6, cex.lab=0.8,cex.axis=0.7));
## also plot for each group
for(i in seq_along(groupVals))
{
	title<-paste(params$main,"(",group,"=",groupVals[i],")",sep="")
	with(subset(dat,get(group)==groupVals[i]),scatter_plot(get(xVar),get(yVar),
	xlim=params$xlim, ylim=params$ylim,
	xlab=params$xlab, ylab=params$ylab,
	main=title, col=plotCols[1+i],
	cex=0.6, cex.lab=0.8,cex.axis=0.7));
}
dev.off()

message("Job is done")

