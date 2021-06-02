#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile="config/config.APA-Scan.yaml" \
    --cores 4 \ # adjust number as needed
    --use-conda \ # or --use-singularity
    --printshellcmds \
    --forcerun generate_configurationfile
