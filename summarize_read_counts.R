#!/usr/bin/env Rscript

stop("Not completed")
res<-try(library("ZZMisc"))
if(class(res) == "try-error") {
    warning("Package 'ZZMisc' isn't installed. Installing it now")
    devtools::install_github("fortune9/ZZMisc")
    library("ZZMisc")
}

load_package("data.table")
load_package("optparse")
load_package("logger")
log_threshold(DEBUG)

# build usage information
desc<-"
Given a bed file, this program summarizes the read counts from input files,
such as bedgraph files, by summing, averaging, ..., of the read counts located
in each region represented in the input bed file.

Options (default value of each option is present in []):
"

epiInfo<-"
Example usage:

%prog <bed-file> <bedgraph-1> [<bedgraph-2> ...]

Author: Zhenguo Zhang
"

exdent<-16

option_list<-list(
                make_option(c("--stat"), action="store", type="character", dest="stat",
                            default="sum", 
                            help=wrap_text("which summary stat to use.",
                                           "Choices: sum, mean.",
                                           "[%default]", exdent=exdent)
                            ),
                make_option(c("-o", "--outfile"), action="store", type="character", dest="outFile",
                            default="",
                            help=wrap_text("The output filename. [%default]")
                  )
)

# parse args
optParser<-OptionParser(
                 usage="usage: %prog [options] <bed-file> <bedgraph-1> [<bedgraph-2> ...]",
                 option_list=option_list,
                 add_help_option=F,
                 description=desc,
                 epilog=epiInfo)

args<-commandArgs(T);

if(length(args) < 2 || any(args %in% c("-h","--help")))
        { print_help(optParser); q("no"); }

log_info("Parsing arguments")
opt<-parse_args2(optParser, args)
posArgs<-opt$args
opt<-opt$options


