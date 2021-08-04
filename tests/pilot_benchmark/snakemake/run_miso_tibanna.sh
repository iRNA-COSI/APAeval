#!/bin/bash
snakemake \
    --configfile="configs/config.miso.yaml" \
    --tibanna \
    --default-remote-prefix="<bucketname>/<subdir>" \
    --cores 4 \
    --use-singularity \
    --printshellcmds 
