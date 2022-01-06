# Command line script to run Roar on a set of BAM files

help <- c("Usage: Roar.R GTF SAMPLE_TABLE BASE_KEY OUTPUT_TSV [--stranded] [--help] [-h]",
          "GTF - Path to Roar annotation GTF file",
          "SAMPLE_TABLE - Path to sample table CSV file used as input to the Roar execution workflow. The CSV should contain 'sample_name', 'bam' and 'condition' columns.",
          "BASE_KEY - Name of 'base' or 'control' condition in the 'condition' column of the sample table",
          "OUTPUT_TSV - Name of/path to output TSV file to write Roar results table",
          "--stranded - (Optional) whether input RNA-seq data were generated with a stranded protocol",
          "--help / -h - (Optional) print this help message and exit")

cl_args <- commandArgs(trailingOnly = TRUE)

if ((length(cl_args) == 0) | ("--help" %in% cl_args) | ("-h" %in% cl_args)) {
  cat(help, sep = "\n")
  stop()
}

library(roar)

gtf_path <- cl_args[1]
sample_tbl_path <- cl_args[2]
base_key <- cl_args[3]
output_tsv <- cl_args[4]

if ("--stranded" %in% cl_args) {
  stranded <- TRUE
} else {
  stranded <- FALSE
}

#1. Extract lists of bams per condition from sample table
sample_tbl <- read.table(sample_tbl_path,
                         header = T,
                         sep = ",",
                         stringsAsFactors = F)


n_cond <- length(unique(sample_tbl$condition))

if ( n_cond != 2) {

  stop(paste("Sample table must only contain 2 distinct conditions -",
             n_cond,
             "were found", ))
}


if (!(base_key %in% sample_tbl$condition)) {

  stop(paste("'base-condition-key' -",
             base_key,
             "- not found in sample table -",
             paste(unique(sample_tbl$condition), sep = ","),
             "were found",
             sep = " ")
  )
}

base <- sample_tbl[sample_tbl$condition == base_key, ]

treat <- sample_tbl[sample_tbl$condition != base_key, ]

# Create lists of paths to bam files for each condition
base_bams <- list(base$bam)
names(base_bams) <- base$sample_name

treat_bams <- list(treat$bam)
names(treat_bams) <- treat$sample_name

# 2. Run Roar

rds <- RoarDatasetFromFiles(treatmentBams = treat_bams,
                            controlBams = base_bams,
                            gtf = gtf_path
)

rds <- countPrePost(rds, stranded = stranded)

rds <- computeRoars(rds)

rds <- computePvals(rds)

results <- totalResults(rds)

# Gene IDs are stored as rownames so add as 1st column
results <- cbind(gene_id = rownames(results), results)

write.table(results,
            output_tsv,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
