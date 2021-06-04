#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.DaPars2.yaml" \
    --cores 4 \
    --use-conda \
    --printshellcmds
# adjust number as needed
# or --use-singularity