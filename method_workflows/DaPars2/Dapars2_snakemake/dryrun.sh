#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.DaPars2.yaml" \
    --cores 4 \
    --use-singularity \
    --printshellcmds \
    --dryrun
