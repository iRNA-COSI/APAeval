#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "config/config.DaPars2.yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.DaPars2.png
