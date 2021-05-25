#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflows/miso.smk" \
    --configfile="configs/config.miso.yaml" \
    --cores 4 \
    --use-conda --conda-frontend conda \
    --printshellcmds 
