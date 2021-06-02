#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflows/miso.smk" \
    --configfile="configs/config.miso.yaml" \
    --cores 4 \
    --use-singularity \
    --singularity-args="--bind ${PWD}/../../test_data"\
    --printshellcmds 