if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("APAlyzer"))) == FALSE ) { stop("[ERROR] Package 'APAlyzer' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("repmis"))) == FALSE ) { stop("[ERROR] Package 'repmis' required! Aborted.") }
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
    "--input_gtf",
    action="store",
    type="character",
    default=FALSE,
    help="Input gtf. Required!",
    metavar="files"
  ),
  make_option(
    "--dir_path",
    action="store",
    type="character",
    default=FALSE,
    help="output path",
    metavar="files"
  ),
  make_option(
    "--out_reference",
    action="store",
    type="character",
    default=FALSE,
    help="Reference genome file compatible with APAlyzer",
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

# Parse command-line arguments
opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] \n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)
pwd <- getwd()
GTFfile <- file.path(pwd, opt$input_gtf)
setwd(opt$dir_path)
pwd1 <- getwd()
file.copy(GTFfile, pwd1)
## build Reference ranges for 3'UTR PASs in mouse	  
PASREFraw=PAS2GEF(GTFfile)
refUTRraw=PASREFraw$refUTRraw
dfIPAraw=PASREFraw$dfIPA
dfLEraw=PASREFraw$dfLE
setwd(pwd)
save(list = c("refUTRraw", "dfIPAraw", "dfLEraw"), file=opt$out_reference)
