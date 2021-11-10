#!/usr/bin/env Rscript

# check the input.csv file
library(rlang)
library(rtracklayer)

args <- commandArgs(trailingOnly = T)
input.df <- read.csv(args[1])

err <- 0

## check header
header <- c("folder_to_BAM","gtf_or_gff","polyA_bed","comparison_table",
            "method","strand","SE_PE","analysis")
if (!all(colnames(input.df) == header)) {
  format_error_bullets(c("x"=paste0("The header must be ", 
                                    paste(sQuote(header), collapse = ","))))
  quit()
}

## check missing values
if (any(is.na(input.df[1,]))) {
  na.cols <- colnames(input.df)[which(is.na(input.df[1,]))]
  format_error_bullets(c("x"=paste0("Argument ", sQuote(na.cols), " empty")))
  err <- err + 1
}

## check the comparison table 
#tb <- read.table(input.df$comparison_table, header = TRUE)

## check some names
if (!grepl("bed", input.df$polyA_bed)) {
  format_error_bullets(c("x"="Please point 'polyA_bed' to a bed file"))
  err <- err + 1
}

if (!input.df$method %in% c("DEXSeq", "limma", "edgeR")) {
  format_error_bullets(c("x"="Please specify 'method' to one of these 3 methods: DEXSeq/limma/edgeR"))
  err <- err + 1
}

if (!input.df$strand %in% c("unstranded", "stranded", "rev_stranded")) {
  format_error_bullets(c("x"="Please specify 'strand' to one of these 3 types: unstranded/stranded/rev_stranded"))
  err <- err + 1
}

if (!input.df$SE_PE %in% c("SE", "PE")) {
  format_error_bullets(c("x"="Please specify 'SE_PE' to either 'SE' or 'PE'"))
  err <- err + 1
}

if (!input.df$analysis %in% c("UTR", "exon")) {
  format_error_bullets(c("x"="Please specify 'analysis' to either 'UTR' or 'exon'"))
  err <- err + 1
}

if (err > 0) {
  quit("Correct above error(s)")
} else {
  if (grepl(".gz", input.df$gtf)) {
    system(paste0("gunzip ", input.df$gtf))
  }
  
  # store a rds file for the polyA atlas bed, otherwise would get an error
  extraColvect <- c(percentage="numeric", numberofprot="integer",
                    tpm2="numeric",encod="character",addinfo="character")
  polyA <- rtracklayer::import(input.df$polyA_bed, 
                               format="bed", extraCols=extraColvect)
  saveRDS(polyA, args[2])
}
