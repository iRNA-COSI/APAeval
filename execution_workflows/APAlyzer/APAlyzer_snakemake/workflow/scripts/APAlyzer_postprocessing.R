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
    "--out_final",
    action="store",
    type="character",
    default=FALSE,
    help="Path for the final output file.",
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
load(opt$out_main)

# Create a dictionary with gene id and p pvalue
out_dict = hash()
for(i in 1:length(out_df)){
    gene_symbol = out_df[i, 1]
    gene_id = gene_dict[[gene_symbol]]
    p_value = out_df[i, 2]
    if(!has.key(gene_id, out_dict) {
        out_dict[[gene_id]] = p_value
    }
    else {
        out_dict[[gene_id]] = min(c(out_dict[[gene_id]], p_value))
    }
}
# Populate a df with gene id and pvalue from dictionary
gene_ids = keys(hash)
pvalues = values(hash)
df = as.data.frame(cbind(gene_ids, pvalues))
df.df = data.frame(lapply(df, as.character), stringsAsFactors=FALSE)

# Writing gene id and pvalue to tsv output file
write.table(mtcars, file = opt$out_final, sep = "\t",
            row.names = FALSE, col.names = NA)

