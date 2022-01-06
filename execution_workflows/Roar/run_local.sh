#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="workflow/config/config.Roar.yaml" \
    --cores 4 \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../../tests/test_data" \
    --printshellcmds \
    
