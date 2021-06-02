#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.APA-Scan.yaml" \
    --cores 4 \
    --use-conda \
    --printshellcmds \
    --dryrun
