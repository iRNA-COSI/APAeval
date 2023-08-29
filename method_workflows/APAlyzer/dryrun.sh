#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.APAlyzer.yaml" \
    --cores 4 \
    --use-singularity \
    --printshellcmds \
    --dryrun