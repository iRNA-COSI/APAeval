#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "config/config.PAQR.yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.PAQR.png
