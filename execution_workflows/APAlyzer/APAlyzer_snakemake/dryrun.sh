#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.[METHOD].yaml" \
    --cores 4 \
    --use-conda \
    --printshellcmds \
    --dryrun
