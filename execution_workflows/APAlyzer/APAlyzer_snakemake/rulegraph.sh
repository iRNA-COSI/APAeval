#!/bin/bash
snakemake \
    --snakefile="workflow/Snakefile" \
    --configfile "config/config.[METHOD].yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.[METHOD].png
