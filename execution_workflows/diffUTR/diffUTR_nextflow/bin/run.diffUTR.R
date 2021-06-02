#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(diffUTR)
  library(SummarizedExperiment)
  library(rtracklayer)
  library(utils)
})

# extract arguments
args <- commandArgs(trailingOnly = T)
input.df <- read.csv(args[1])
bed.rds <- args[2]
out.file <- args[3]

bamfiles <- list.files(input.df$folder_to_BAM, pattern = "*.bam$", full.names=TRUE)
gtf <- gsub(".gz$", "", input.df$gtf)
comp.df <- read.table(input.df$comparison_table)
# 0-unstranded; 1-stranded; 2-reversely stranded
strand.info <- ifelse(input.df$strand=="unstranded", 0, 
                      ifelse(input.df$strand=="stranded", 1, 2))
is.PE <- ifelse(input.df$SE_PE=="PE", TRUE, FALSE)
if (input.df$method == "DEXSeq") {
  diff <- utils::getFromNamespace("DEXSeqWrapper", "diffUTR")
} else if (input.df$method == "limma") {
  diff <- utils::getFromNamespace("diffSpliceWrapper", "diffUTR")
} else {
  diff <- utils::getFromNamespace("diffSpliceDGEWrapper", "diffUTR")
}


# run diffUTR
bins <- prepareBins(g=gtf, APA=bed.rds)
rse <- countFeatures(bamfiles, bins, strandSpecific=strand.info, 
                     nthreads=12, isPairedEnd=is.PE)

colData(rse) <- DataFrame(comp.df)
rse <- diff(rse, design = ~ condition)
utr <- geneLevelStats(rse, includeTypes="UTR", returnSE=FALSE)

# store output file
utr2 <- as.data.frame(utr)
gene.id <- rownames(utr2)
utr2 <- cbind.data.frame(geneID=gene.id, utr2)
write.table(utr2, file=out.file, quote=FALSE, row.names=FALSE, 
            col.names = TRUE, sep="\t")
