if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("APAlyzer"))) == FALSE ) { stop("[ERROR] Package 'APAlyzer' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("repmis"))) == FALSE ) { stop("[ERROR] Package 'repmis' required! Aborted.") }
# if ( suppressWarnings(suppressPackageStartupMessages(require("diffloop"))) == FALSE ) { stop("[ERROR] Package 'diffloop' required! Aborted.") }

#######################
###  PARSE OPTIONS  ###
#######################

# Get script name
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Build description message
description <- "Build reference genome\n"
version <- "Version: 1.0.0 (May 2021)"
requirements <- "Requires: optparse, APAlyzer"
msg <- paste(description, version, requirements, sep="\n")

# Define list of arguments
option_list <- list(
  make_option(
    "--reference_genome",
    action="store",
    type="character",
    default=FALSE,
    help="Input gtf. Required!",
    metavar="files"
  ),
  make_option(
    "--sample_name",
    action="store",
    type="character",
    default=FALSE,
    help="Input samples table. Required!",
    metavar="files"
  ),
  make_option(
    "--sample_path",
    action="store",
    type="character",
    default=FALSE,
    help="Input samples table. Required!",
    metavar="files"
  ),
  make_option(
    "--read_orientation",
    action="store",
    type="character",
    default=FALSE,
    help="Input samples table. Required!",
    metavar="files"
  ),
  # StrandType="forward-reverse"    ## "forward-reverse",  or "reverse-forward" or "NONE"	
  make_option(
    "--out_quantification",
    action="store",
    type="character",
    default=FALSE,
    help="output quantification file",
    metavar="files"
  ), 
  make_option(
    c("-h", "--help"),
    action="store_true",
    default=FALSE,
    help="Show this information and die."
  ),
  make_option(
    c("-v", "--verbose"),
    action="store_true",
    default=FALSE,
    help="Print log messages to STDOUT."
  )
)

opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] \n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

flsall<-opt$sample_path 
names(flsall)<-opt$sample_name

read_orientation <-opt$read_orientation

load(opt$reference_genome)

PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
UTRdbraw=PASREF$UTRdbraw
dfIPA=PASREF$dfIPA
dfLE=PASREF$dfLE
print('UTRdbraw')
print(head(UTRdbraw))
print(class(UTRdbraw))
print('dfIPA')
print(head(dfIPA))
print(class(dfIPA))
print('dfLE')
print(head(dfLE))
print(class(dfLE))
IPA_OUT=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype=read_orientation, nts=1)

print('ok5')
write.table(IPA_OUT, file = opt$out_quantification, row.names=FALSE)
print('ok6')

