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
    "--dir_path",
    action="store",
    type="character",
    default=FALSE,
    help="output path",
    metavar="files"
  ),
  make_option(
    "--out_main",
    action="store",
    type="character",
    default=FALSE,
    help="Preprocessing variables needed for APAlyzer",
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

# Load variables from preprocessing step
load(opt$out_preprocessing)

#----------------------Build the PAS reference regions----------------------
PASREFraw=PAS2GEF(GTFfile, AnnoMethod="V2")
refUTRraw=PASREFraw$refUTRraw
dfIPAraw=PASREFraw$dfIPA
dfLEraw=PASREFraw$dfLE
PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
# reference region used for 3'UTR APA analysis
UTRdbraw=PASREF$UTRdbraw
# dfIPA and dfLE are needed in intronic APA analysis
dfIPA=PASREF$dfIPA
dfLE=PASREF$dfLE

#----------Calculation of relative expression of 3'UTR APA and IPA----------
UTR_APA_OUT=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="invert")
IPA_OUT=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="invert", nts=1)

#------------------Significantly regulated APA in 3’UTRs--------------------
############# 3utr APA #################
sampleTable = data.frame(samplename =
                           names(flsall),
                         condition = condition_counts)

out_3UTRAPA=APAdiff(sampleTable,UTR_APA_OUT,
                     conKET=unique_conditions[1],
                     trtKEY=unique_conditions[2],
                     PAS='3UTR',
                     CUTreads=5)

#-----------------Significantly regulated APA in Intron--------------------
############# IPA #################
out_IPA=APAdiff(sampleTable,
                 IPA_OUT,
                 conKET=unique_conditions[1],
                 trtKEY=unique_conditions[2],
                 PAS='IPA',
                 CUTreads=5)

# merge the output dataframes for postprocessing
out_3UTRAPA = out_3UTRAPA[,c("gene_symbol", "pvalue")]
out_IPA = out_IPA[,c("gene_symbol", "pvalue")]
out_df = rbind(out_3UTRAPA, out_IPA)

# Save the variables needed for APAlyzer_main
setwd(pwd)
gene_dict = gene_dict
save(list = c("gene_dict", "out_df"), file=opt$out_main)
