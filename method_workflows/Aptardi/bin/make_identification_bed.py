#!/usr/bin/env python

import sys
import pyranges as pr
import pandas as pd

input_gtf_path = sys.argv[1]
output_bed_path = sys.argv[2]

if output_bed_path == 'not_specified':
    output_bed_path = input_gtf_path.split(".")[0]+"_identification.bed"

gtf = pr.read_gtf(input_gtf_path)
# Source is 2nd col in GTF, Feature is 3rd col
gtf = gtf.subset(lambda df: (df["Source"] == "aptardi") & (df["Feature"] == "transcript"))

if len(gtf) == 0:
    raise Exception("No Aptardi predicted transcripts found in provided GTF file")

# Strand aware, reports 3'most coordinate for each interval
# PyRanges coordinates follow BED conventions
three_ends = gtf.three_end()
# Some reformatting to get a BED6 file
three_ends = three_ends[["Strand"]]
# Dummy score column
three_ends.Score = "."

# Drop duplicate polyA sites if present
three_ends = three_ends.drop_duplicate_positions()

# Making up a Name column, can be whatever you like though
# <chr>:<start>:<end>
three_ends = three_ends.assign("Name",
                               lambda df: df[["Chromosome", "Start", "End"]].apply(lambda row: ":".join(row.values.astype(str)), axis=1))

three_ends.to_bed(output_bed_path)
