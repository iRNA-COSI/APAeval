#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "config/config.Roar.yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.Roar.png
