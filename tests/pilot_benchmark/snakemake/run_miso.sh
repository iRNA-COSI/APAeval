#!/bin/bash
snakemake \
    --rerun-incomplete \
    --snakefile="Snakefile" \
    --configfile="configs/config.miso.yaml" \
    --cores 4 \
    --use-conda --conda-frontend conda \
    --printshellcmds --dryrun
