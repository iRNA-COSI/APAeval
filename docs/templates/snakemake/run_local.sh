#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.[METHOD].yaml" \
    --cores 4 \ # adjust number as needed
    --use-conda \ # or --use-singularity
    --printshellcmds
