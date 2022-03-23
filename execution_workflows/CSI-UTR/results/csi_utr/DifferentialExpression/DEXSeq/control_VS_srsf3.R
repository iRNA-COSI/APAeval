#############################################################
## CSI-UTR srsf3 vs control ##
## SCRIPT AUTOMATICALLY GENERATED                          ##
#############################################################
if (!require("DESeq2")) { 
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2", dependencies=TRUE)  ## INSTALL DESeq2
    if(!require("DESeq2")) stop("DESeq2 package not found")
}

if (!require("DEXSeq")) { 
    source("http://bioconductor.org/biocLite.R")
    biocLite("DEXSeq", dependencies=TRUE)  ## INSTALL DEXSeq
    if(!require("DEXSeq")) stop("DEXSeq package not found")
}

library("DESeq2")
library("DEXSeq")

setwd("/home/fzhuang/apaeval/APAeval/execution_workflows/CSI-UTR/./results/csi_utr//DifferentialExpression/DEXSeq/control_VS_srsf3/")
fileArr = c("../../../DEXSeq/siSrsf3_R1_2genes_Chr_prefix.DEXSeq.counts", "../../../DEXSeq/siControl_R1_2genes_Chr_prefix.DEXSeq.counts")
gffFile = "../../../DEXSeq/siSrsf3_R1_2genes_Chr_prefix.DEXSeq.gff"
sampleTable = data.frame(
   row.names=c("siSrsf3_R1_2genes_Chr_prefix", "siControl_R1_2genes_Chr_prefix"),
   condition=c("A_srsf3", "B_control"),
   libType=c("single-end", "single-end")
)

diffExp = DEXSeqDataSetFromHTSeq (fileArr, sampleData=sampleTable, 
                                  design=~sample+exon+condition:exon,
                                  flattenedfile=gffFile)

#########################################
## MAY WANT TO UNCOMMENT FOR DEBUGGING ##
#########################################

# colData(diffExp)
# head(counts(diffExp), 5)
# head(featureCounts(diffExp), 5)
# sampleAnnotation(diffExp)

diffExp = estimateSizeFactors(diffExp)  ## NORMALIZE DATA 
diffExp = estimateDispersions(diffExp)  ## ESTIMATE DISPERSION 
# plotDispEst(diffExp) ## PRINT DISPERSION ESTIMATES 
diffExp = testForDEU(diffExp) ## TEST FOR SIGNIFICANT EVENTS
diffExp = estimateExonFoldChanges(diffExp) ## ADD FOLD CHANGE INFORMATION
diffExpRegions = DEXSeqResults(diffExp)

sigP <- diffExpRegions[!(is.na(diffExpRegions$pvalue)), ]
defExp <- sigP
sigP <- sigP[(sigP$pvalue < 0.05), ]
sigFDR <- sigP[!(is.na(sigP$padj)), ]
sigFDR <- sigFDR[(sigFDR$padj < 0.10), ]
defExp <- defExp[order(defExp$pvalue), ]
sigFDR <- sigFDR[order(sigFDR$padj), ]

# defExp
# sigFDR

write.table(defExp, file="allCSIs_control_VS_srsf3_SORTED.txt", sep="\t")
write.table(diffExpRegions, file="allCSIs_control_VS_srsf3_BYREGION.txt", sep="\t")
write.table(sigFDR, file="sigCSIs_control_VS_srsf3.txt", sep="\t")


# PLOT RESULTS FOR CAMK4
# plotDEXSeq(dxr1, "ENSMUSG00000038128", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)

