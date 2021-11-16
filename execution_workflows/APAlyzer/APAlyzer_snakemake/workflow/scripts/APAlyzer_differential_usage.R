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
    "--sample_table",
    action="store",
    type="character",
    default=FALSE,
    help="Input table with samples. Required!",
    metavar="files"
  ),
  make_option(
    "--outfile",
    action="store",
    type="character",
    default=FALSE,
    help="output differential usage file",
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



opt_parser <- OptionParser(usage=paste("Usage:", script, "[OPTIONS] --input <paths/to/input/tables>\n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

sampleTable<-read.csv(opt$sample_table, header = TRUE)
sampleTable<-sampleTable[,c("sample", "condition")]

############# 3utr APA #################
# sampleTable = data.frame(samplename = 
# c('Heart_rep1',
# 'Heart_rep2',
# 'Heart_rep3',
# 'Heart_rep4',
# 'Testis_rep1',
# 'Testis_rep2',
# 'Testis_rep3',
# 'Testis_rep4'),
# condition = c(rep("Heart",4),
# rep("Testis",4)))


						
test_3UTRAPA=APAdiff(sampleTable,opt$UTR_APA_OUT, 
                        conKET='control',
                        trtKEY='treatment',
                        PAS='3UTR',
                        CUTreads=5)
						
table(test_3UTRAPA$APAreg)
write.csv(UTR_APA_OUT, opt$outfile, row.names=FALSE)				
APAVolcano(test_3UTRAPA, PAS='3UTR', Pcol = "pvalue", plot_title='3UTR APA')
APABox(test_3UTRAPA, xlab = "APAreg", ylab = "RED", plot_title = NULL)