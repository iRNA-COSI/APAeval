
### Libraries

library(roar)
library(rtracklayer)
library(org.Mm.eg.db)

# Input files
args = commandArgs(trailingOnly=TRUE)

# Input gtf
gtf <- args[1]

# Input bams
bam1 <- paste(path=args[2],list.files(path=args[2],pattern="*.bam"),sep="/")
bam2 <- paste(path=args[3],list.files(path=args[3],pattern="*.bam"),sep="/")
bams <- c(bam1, bam2)
# Output prefix
out.prefix <- args[4]






# Print out inputs
print(">>>>> INPUT ARGS ARE:")
print(paste0("   >> GTF:  ",args[1]))
print(paste0("   >> BAM1: ",bam1))
print(paste0("   >> BAM2: ",bam2))




### Create RoarData Object

# Generate RoarDataSet
print(">>>>> GENERATING ROAR DATA SET")
rds <- RoarDatasetMultipleAPAFromFiles(bams[1:2], bams[3:4], gtf)




### Compute Roar and Associated Statistics/Metrics

# UTR Counts
print(">>>>> COUNTING 3UTR READS")
rds <- countPrePost(rds, FALSE)

# Compute Roars
print(">>>>> COMPUTING ROAR")
rds <- computeRoars(rds)

# Compute p-values, using Fisher's Exact test, considering pairing of samples and combining p-values
# for those of the same level
print(">>>>> COMPUTING PAIRED PVALUES")
rds <- computePairedPvals(rds, c(1,2), c(1,2))




### Get Output

# Get results df
print(">>>>> GETTING OUTPUT")
results.paired <- totalResults(rds)

# Add columns for UTR ids & gene IDS
results.paired["3utrid"] <- gsub(".*_","",rownames(results.paired))
results.paired["ENTREZID"] <- gsub("_.*","",rownames(results.paired))




### Convert Gene IDs

# Get ENTREZID
print(">>>>> CONVERTING IDS")
geneids <- gsub("_.*","",rownames(results.paired))

# Map to ENSEMBL gene Id
enz2ens <- select(org.Mm.eg.db, # The annotation library object (contains the info needed to convert ids)
                   keys=geneids,
                   columns=c("ENTREZID","ENSEMBL"),
                   keytype="ENTREZID")

# Get gtf as GRanges and format to PolyA site
gtf.granges <- import(gtf)
gtf.granges$apa <- gsub("_.*","",gtf.granges$apa)
gtf.granges <- data.frame(gtf.granges,stringsAsFactors=FALSE)
colnames(gtf.granges)[[10]] <- "3utrid"

# Merge it all together
outdf <- merge(enz2ens,results.paired,by="ENTREZID")
outdf <- merge(gtf.granges,outdf,by="3utrid")




### Format output

print(">>>>> SAVING RESULTS")
# Output df - format 1 (BED) & format 3
outdf01 <- outdf[,c("seqnames","start","end","source","3utrid","strand")]; outdf01["source"] <- "."
outdf03 <- outdf[,c("ENSEMBL","pval")]; outdf03 <- outdf03[!duplicated(outdf03),]

# save files
write.table(outdf01,paste0(out.prefix,"_01.bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(outdf03,paste0(out.prefix,"_03.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
