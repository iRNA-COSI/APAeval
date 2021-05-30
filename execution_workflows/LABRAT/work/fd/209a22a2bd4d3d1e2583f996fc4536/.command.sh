#!/bin/bash -euo pipefail
LABRAT.py --mode makeTFfasta --gff test.gff3 --genomefasta mm39.fa --lasttwoexons --librarytype RNAseq
