#!/usr/bin/env Rscript

## functions
# get the script path
self<-function()
{
        args<-commandArgs(F)
        path<-sub("^--file=",'',grep("^--file=", args, value=T))
                normalizePath(path)
}

read_data<-function(f, ...) {
    dat<-fread(f, ...)
    return(dat)
}

write_data<-function(dat, ...) {
    fwrite(dat, ...)
}

# for each region, get a string of all overlapped features
get_ids<-function(rowIds, dat, feat) {
	if(is.na(rowIds[1])) { return("") }
	ids<-unique(unlist(dat[rowIds, feat, with=F]))
	return(paste(ids, collapse=";"))
}

library("optparse")
library("data.table")

desc="
This program adds overlapped features to query
regions.
"
epiInfo="
Example usage:
%prog <query-region-file> <feature-file>
Author: Zhenguo Zhang
Updated: Tue Aug 31 05:37:02 UTC 2021
"

option_list<-list(
           make_option("--pos-col",
                       type="integer", action="store", default=NULL,
                       dest="posCol",
                       help="the column where genomic positions are stored. If this is given, the start and end positions are infered from it"),
           make_option("--feat-name",
                       type="character", action="store", default="feature",
                       dest="featName",
                       help="the name of column for the added feature in output [%default]"),
           make_option("--out-file",
                       type="character", action="store", default="",
                       dest="outFile",
                       help="the output filename [console]")
           )

optParser<-OptionParser(
                        usage="usage: %prog [options] <query-region-file> <feature-file>",
                        option_list=option_list,
                        description=desc,
                        epilog=epiInfo
                        )

args<-commandArgs(T);

if(length(args) < 2 || any(args %in% c("-h","--help")))
        { print_help(optParser); q("no"); }

opt<-parse_args2(optParser, args)
posArgs<-opt$args
opt<-opt$options

if(length(posArgs) < 2 ) {
    stop("Both <query-region-file> and <feature-file> are required")
}

posCol<-opt$posCol
outFile<-opt$outFile
featName<-opt$featName

message("Reading data")
regions<-read_data(posArgs[1])
feats<-read_data(posArgs[2], header=F)

# set key for features
setkey(feats, V1, V2, V3)
feats<-feats[V2 <= V3,] # remove start > end

# overlap
if(is.null(posCol)) { # not single position
    regionCols<-colnames(regions)[1:3]
} else {
    # convert 1-based position to 0-based start and end, and cover CG
    # dinucleotide
    message("Coverting 1-based positions to half-open regions")
    regions[, startZZ:=regions[[posCol]]-1]
    regions[, endZZ:=regions[[posCol]]+1]
    regionCols<-c(colnames(regions)[1], "startZZ", "endZZ")
}

message("Overlapping regions")
op<-foverlaps(regions, feats, 
              by.x=regionCols, by.y=c("V1","V2","V3"), 
              type="any", mult="all", nomatch="NA", which=T)
## collapse all overlapped features for each region
featIds<-op[, get_ids(yid, feats, "V4"), by=xid]
#save.image("tmp.RDa"); q("no")
regions[, (featName) := featIds$V1]
if("startZZ" %in% colnames(regions)) {
    regions$startZZ<-NULL
    regions$endZZ<-NULL
}

message("Writing results")
write_data(regions, file=outFile, sep="\t", na='""')

message("Job is done")

