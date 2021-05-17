#!/usr/bin/bash
snakemake \
    --configfile "config/config.[METHOD].yaml" \
    --rulegraph -np | dot -Tpng > rulegraph.[METHOD].png
