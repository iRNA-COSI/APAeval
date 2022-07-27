#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflow/Snakefile" \
    --configfile="workflow/config/config.Roar.yaml" \
    --cores 4 \
    --use-singularity \
    --printshellcmds \
    --dryrun
