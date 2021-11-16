#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "config/config.APAlyzer.yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.APAlyzer.png
