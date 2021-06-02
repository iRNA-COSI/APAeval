#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "config/config.APA-Scan.yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.APA-Scan.png
