#!/bin/bash
snakemake \
    --configfile="configs/config.miso.yaml" \
    --cores 4 \
    --use-conda --conda-frontend conda \
    --printshellcmds --dryrun
