#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.DaPars2.yaml" \
    --cores 4 \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../../../tests/test_data" \
    --printshellcmds
# adjust number of cores as needed
