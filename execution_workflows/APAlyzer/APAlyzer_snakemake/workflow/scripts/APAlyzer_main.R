#------------------------Main---------------------------
# Run APAlyzer using variables from preprocessing step

# load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("APAlyzer"))) == FALSE ) { stop("[ERROR] Package 'APAlyzer' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("GenomicFeatures"))) == FALSE ) { stop("[ERROR] Package 'repmis' required! Aborted.") }

#########################
###  PARSE ARGUMENTS  ###
#########################

# Get script name
script = sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

# Build description message
description = "Build reference genome\n"
version = "Version: 1.0.0 (May 2021)"
requirements = "Requires: optparse, APAlyzer, GenomicFeatures"
msg = paste(description, version, requirements, sep="\n")

# Define list of arguments
option_list = list(
  make_option(
    "--dir_path",
    action = "store",
    type = "character",
    default = FALSE,
    help = "output path",
    metavar = "files"
  ),
  make_option(
    "--in_main",
    action = "store",
    type = "character",
    default = FALSE,
    help = "Input from preprocessing step to run APAlyzer",
    metavar = "files"
  ),
  make_option(
    "--out_main",
    action = "store",
    type = "character",
    default = FALSE,
    help = "APAlyzer output variables to postprocess",
    metavar = "files"
  ), 
  make_option(
    c("-h", "--help"),
    action = "store_true",
    default = FALSE,
    help = "Show this information and die."
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print log messages to STDOUT."
  )
)

# Parse command-line arguments
opt_parser = OptionParser(
                usage = paste("Usage:", script, "[OPTIONS] \n", sep = " "),
                option_list = option_list,
                add_help_option = FALSE,
                description = msg)
opt = parse_args(opt_parser)

######################
###  RUN APALYZER  ###
######################

# Load variables from preprocessing step
load(opt$in_main)

#----------------------Build the PAS reference regions----------------------
PASREFraw = PAS2GEF(GTFfile, AnnoMethod="V2")
refUTRraw = PASREFraw$refUTRraw
dfIPAraw = PASREFraw$dfIPA
dfLEraw = PASREFraw$dfLE
PASREF = REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
# reference region used for 3'UTR APA analysis
UTRdbraw = PASREF$UTRdbraw
# dfIPA and dfLE are needed in intronic APA analysis
dfIPA = PASREF$dfIPA
dfLE = PASREF$dfLE

#----------Calculation of relative expression of 3'UTR APA and IPA----------
# calculate relative expression of 3' UTR APA
UTR_APA_OUT = PASEXP_3UTR(UTRdbraw, flsall, Strandtype = "NONE")

# ensure that coordinates are numeric
dfIPA$Pos = as.numeric(as.character(dfIPA$Pos))
dfIPA$upstreamSS = as.numeric(as.character(dfIPA$upstreamSS))
dfIPA$downstreamSS = as.numeric(as.character(dfIPA$downstreamSS))
dfLE$LEstart = as.numeric(as.character(dfLE$LEstart))
dfLE$TES = as.numeric(as.character(dfLE$TES))

# calculate relative expression of IPA
IPA_OUT = PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype = "NONE", nts = 1)

#---------------Significantly regulated APA in 3’UTRs and introns-----------------
sampleTable = data.frame(samplename =
                           names(flsall),
                         condition = condition_counts)

if(nrow(UTR_APA_OUT) > 1) {
    out_3UTRAPA = APAdiff(sampleTable,UTR_APA_OUT,
                        conKET = unique_conditions[1],
                        trtKEY = unique_conditions[2],
                        PAS = '3UTR',
                        CUTreads = 5)
    out_3UTRAPA = out_3UTRAPA[,c("gene_symbol", "pvalue")]
} else {
    out_3UTRAPA = NULL	
}

if(nrow(IPA_OUT) > 1) {
    out_IPA = APAdiff(sampleTable,
                    IPA_OUT,
                    conKET = unique_conditions[1],
                    trtKEY = unique_conditions[2],
                    PAS = 'IPA',
                    CUTreads = 5)
    out_IPA = out_IPA[,c("gene_symbol", "pvalue")]
} else {
   out_IPA = NULL	
}

############################
###  SAVE FINAL OUTPUTS  ###
############################

# merge apalyzer output dataframes for postprocessing
out_df = rbind(out_3UTRAPA, out_IPA)

# Save the variables needed for APAlyzer_postprocessing
pwd = getwd()
setwd(pwd)
gene_dict = gene_dict
save(list = c("gene_dict", "out_df"), file = opt$out_main)
