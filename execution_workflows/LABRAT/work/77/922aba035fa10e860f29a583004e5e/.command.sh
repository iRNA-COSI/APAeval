#!/bin/bash -euo pipefail
conda run -n labrat LABRAT.py --mode runSalmon --librarytype RNAseq --txfasta TFseqs.fasta --reads1 /home/wanyk/Downloads/APAeval_local/test_data/siSrsf3_R1_100genes.fastq.gz --samplename test --threads 8
