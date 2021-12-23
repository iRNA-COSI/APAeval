#------------------------Preprocessing------------------------
# 1. Check sample file column names
# 2. Check that number of conditions is exactly 2
# 3. Get list of sample names and absolute paths to bam file
# 4. Get a dictionary of gene symbol to gene id from gtf file

if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("stringr"))) == FALSE ) { stop("[ERROR] Package 'stringr' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("hash"))) == FALSE ) { stop("[ERROR] Package 'hash' required! Aborted.") }

# Build description message
description <- "Build reference genome\n"
version <- "Version: 1.0.0 (May 2021)"
requirements <- "Requires: optparse, stringr, hash"
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
    "--sample_file_path",
    action="store",
    type="character",
    default=FALSE,
    help="Sample file path. Required!",
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
    "--out_preprocessing",
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
opt_parser = OptionParser(usage=paste("Usage:", script, "[OPTIONS] \n", sep=" "), option_list = option_list, add_help_option=FALSE, description=msg)
opt = parse_args(opt_parser)
pwd = getwd()
GTFfile = file.path(pwd, opt$input_gtf)
setwd(opt$dir_path)
pwd1 = getwd()
file.copy(GTFfile, pwd1)

# Read the sample file
df = read.csv(file=opt$sample_file_path)

# Check column names of the sample file
colnames = c("sample", "bam", "condition", "orientation")
assert_colnames(df, colnames, only_colnames = TRUE, quiet = FALSE)

# Check that there are exactly two conditions
conditions = df[, "condition"]
unique_conditions = conditions[!duplicated(conditions)]
if(length(unique_conditions != 2)) {
    print(paste0("Number of conditions in sample file should be exactly 2, got ", length(conditions)), quote=FALSE)
    exit()
}

# Rearrange rows in sample file to group by condition
df = df[order(df$condition),]

# Get list of sample names
flsall = df[, "bam"]
names(flsall) = df[, "sample"]

# Get list of conditions and counts
condition_counts = c()
for(condition in conditions){
    condition_counts = c(condition_counts, rep(condition, sum(conditions==condition))
}

# Get a dictionary of gene symbol to gene id from gtf file
# initialize a dictionary
gene_dict = hash()
for(row in df[,9]) {
  for(str in str_split(row, ';')[[1]]){
    str = trimws(str)
    gene_symbol = ""
    # get gene id
    key_value_pair = str_split(str, " ")[[1]]
    key = key_value_pair[1]
    value = key_value_pair[2]
    if(key == "gene_id"){
      gene_id = value
    }
    # get gene symbol
    if(key== "gene_name"){
      gene_symbol = value
    }
    if(gene_symbol != "") {
      gene_dict[[gene_symbol]] = gene_id
    }
  }
}

# Save the variables needed for APAlyzer_main
setwd(pwd)
save(list = c("flsall", "gene_dict", "GTFfile", "unique_conditions", "condition_counts"), file=opt$out_preprocessing)