#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.APAlyzer.yaml" \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../../../tests/test_data" \
    --cores 4 \
    --printshellcmds
# adjust number of cores as needed
