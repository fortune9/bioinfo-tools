#!/bin/env Rscript
# functions
region_overlap<-function(v) # input is (start1, end1, start2, end2)
{
	if(is.na(v[3])) {return(NA)}
	s.max<-ifelse(v[1] < v[3], v[3], v[1])
	e.min<-ifelse(v[2] < v[4], v[2], v[4])
	return(e.min-s.max+1)
}

region_overlap1<-function(s1,e1,s2,e2) # input is (start1, end1, start2, end2)
{
	if(is.na(s2)) {return(as.numeric(NA))}
	s.max<-ifelse(s1 < s2, s2, s1)
	e.min<-ifelse(e1 < e2, e1, e2)
	return(e.min-s.max+1)
}

library(data.table)

#eFile<-"sample1.read_loc.tsv.gz"
#oFile<-"in755_1.read_loc.tsv.gz"
inputs<-commandArgs(T)
if(length(inputs) < 3)
{
	stop("Usage: determine_align_accuracy.R expected-file.gz observed-file.gz outfile-base")
}

eFile<-inputs[1];
oFile<-inputs[2];
outbase<-inputs[3];

message("Reading data\n")
# read into data for each chromosome and sum up together
results<-data.frame()
localChr<-"chr1"
expected<-fread(paste("zcat ", eFile, " | sed -ne '2,$p' | sed -e 's/[:-]/\\t/g' | gawk '$1==\"", localChr, "\" && NF==4'", sep=""))
observed<-fread(paste("zcat ", oFile, " | gawk '$2==\"", localChr, "\"'", sep=""), head=F)
setnames(expected, c("chr.e","start.e","end.e", "id"))
setnames(observed, c("id", "chr.o","start.o","end.o", "mapq"))
tmp1<-merge(expected, observed, by="id", all.x=T)

# get statistics
message("Calculating statistics\n")
num.total<-tmp1[,.N]
num.mapped<-tmp1[!is.na(start.o),.N]
tmp1[,diff.start:=start.o-start.e]
tmp1[,diff.end:=end.o-end.e]
#tmp2<-apply(tmp1[,.(start.e,end.e,start.o,end.o)], 1, region_overlap)
tmp2<-tmp1[,region_overlap1(start.e,end.e,start.o,end.o), by=1:nrow(tmp1)]
#tmp1[, overlap:=tmp2]
tmp1[, overlap:=tmp2[,V1]]
breaks<-c(0,1,2,5,10,100,1e9)
labels<-c("<=1","<=2","<=5","<=10","<=100",">100")
diffCnts<-tmp1[, table(cut(abs(diff.start)+abs(diff.end), breaks, include.lowest=T, right=T, labels=labels))]
num.correct<-diffCnts[1]
breaks<-c(-1e9,0,0.9,0.99,1.1)
labels<-c("<0","0~0.9","0.9~0.99",">0.99")
overlapFrac<-tmp1[, table(cut(overlap/(end.e-start.e), breaks, include.lowest=T, right=T, labels=labels))]

# make plot
message("Making plots\n")
outfig<-paste("outbase","pdf", sep=".")
pdf(file=outfig, width=7, height=4)
layout(matrix(1:2,nr=1))
barplot(diffCnts, las=3)
barplot(overlapFrac, las=3)
dev.off()

cat(sprintf("
 Total number of input reads: %d,
Total number of mapped reads: %d,
      Correctly mapped reads: %d
", num.total, num.mapped, num.correct
))

message("Job done\n")

