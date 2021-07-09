#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# get input table

print(args)
df = read.csv(args[1], header=TRUE, sep =",")
path = args[2]

df = df[,c("bam_prefix", "condition", "sample")]

colnames(df) = NULL

write.table(df, paste0(path,"/sample_table.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
