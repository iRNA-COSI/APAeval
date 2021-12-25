#------------------------Postprocessing------------------------
# Takes APAlyzer output and generates differential challenge
# output tsv file

# load libraries
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("[ERROR] Package 'optparse' required! Aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("hash"))) == FALSE ) { stop("[ERROR] Package 'hash' required! Aborted.") }

#########################
###  PARSE ARGUMENTS  ###
#########################

# Get script name
script = sub("--file=", "", basename(commandArgs(trailingOnly = FALSE)[4]))

# Build description message
description = "Build reference genome\n"
version = "Version: 1.0.0 (May 2021)"
requirements = "Requires: optparse, APAlyzer"
msg = paste(description, version, requirements, sep = "\n")

# Define list of arguments
option_list = list(
  make_option(
    "--in_postprocessing",
    action = "store",
    type = "character",
    default = FALSE,
    help = "APAlyzer output variables to preprocess.",
    metavar = "files"
  ),
  make_option(
    "--out_postprocessing",
    action = "store",
    type = "character",
    default = FALSE,
    help = "Path for the final output file.",
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
                usage=paste("Usage:", script, "[OPTIONS] \n", sep=" "),
                option_list = option_list,
                add_help_option = FALSE,
                description = msg)
opt = parse_args(opt_parser)

#########################
###  GENERATE OUTPUT  ###
#########################

# Load variables from preprocessing step
load(opt$in_postprocessing)

# Create a dictionary with gene id and p pvalue
out_dict = hash()
for(i in 1:length(out_df)) {
    gene_symbol = out_df[i, 1]
    gene_id = gene_dict[[gene_symbol]]
    p_value = out_df[i, 2]
    if(!has.key(gene_id, out_dict)) {
        out_dict[[gene_id]] = p_value
    } else {
        out_dict[[gene_id]] = min(c(out_dict[[gene_id]], p_value))
    }
}
# Populate a df with gene id and pvalue from dictionary
gene_ids = keys(out_dict)
pvalues = values(out_dict)
df = as.data.frame(cbind(gene_ids, pvalues))
df.df = data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

# Writing gene id and pvalue to tsv output file
write.table(df, file = opt$out_postprocessing, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

