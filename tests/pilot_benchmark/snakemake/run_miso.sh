#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="workflows/miso.smk" \
    --configfile="configs/config.miso.yaml" \
    --cores 10 \
    --use-conda \
    --printshellcmds --dryrun