#!/bin/bash

# using BEDTools to intersect the file containing gene names and coordinates with the PAS
# chromosome naming convention must match between the two files

bedtools intersect -wo -s -a SRX351950.clusters.2.0.GRCh38.96.bed -b only_genes_hsap.bed > s2.bed
