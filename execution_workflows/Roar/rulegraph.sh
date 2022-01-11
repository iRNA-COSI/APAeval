#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "workflow/config/config.Roar.yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.Roar.png
