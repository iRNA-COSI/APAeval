#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.PAQR.yaml" \
    --cores 4 \
    --use-singularity \
    --singularity-args "--bind $PWD/../../../" \
    --printshellcmds
