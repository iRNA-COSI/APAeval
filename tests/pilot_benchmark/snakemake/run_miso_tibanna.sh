#!/bin/bash
snakemake \
    --configfile="configs/config.miso.yaml" \
    --tibanna \
    --default-remote-prefix="apaevaltestbucket/test" \
    --cores 4 \
    --use-singularity \
    --printshellcmds 
