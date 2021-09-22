#!/bin/bash
snakemake --unlock;
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.APAlyzer.yaml" \
    --use-conda \
    --conda-frontend conda \
    --cores 4 \
    --printshellcmds 
    
# --use-singularity \
# --singularity-args "--bind ../../../" \
