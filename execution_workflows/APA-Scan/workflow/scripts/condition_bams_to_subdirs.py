#!/usr/bin/env python3


import pandas as pd
from itertools import starmap
from fnmatch import fnmatch
import os
import sys
import argparse

"""
APA-Scan requires BAM files for each condition to be in separate directories
This feels a little restrictive at the pipeline configuration step
This script takes as input:
- sample table for the pipeline
- output directory path

And outputs:
- Directories under output_directory_path, 1 for each condition
- All BAM files for samples of the same condition under their respective
"""


def main(sample_tbl, outdir):
    """
    """

    sample_tbl = pd.read_csv(sample_tbl)

    # check column headers
    assert "sample,bam,condition,data_type".split(",") == sample_tbl.columns.tolist(), "Invalid column headers in input sample table, must be - 'sample,bam,condition,data_type'"

    # check that only two input conditions in sample table (max for APA-Scan)
    assert len(set(sample_tbl["condition"])) == 2, f"Sample table must only contain two distinct conditions under 'condition' in sample_tbl - {len(set(sample_tbl['condition']))} distinct conditions found"

    #1. Make subdirectories under outdir corresponding to each condition

    conditions = list(set(sample_tbl["condition"]))

    for cond in conditions:
        subdir = os.path.join(outdir, cond)

        if os.path.isdir(subdir):
            raise Exception(f"{subdir} already exists - exiting to avoid incorporation of .bam files not specified in sample table")

        os.makedirs(subdir)

    # generate list of tuples of (path_to_bam, out_soft_link_path) for each sample
    bam_in_out_paths = [(os.path.abspath(bam),
                        os.path.join(outdir, condition, sample + ".bam")
                        ) for sample, bam, condition in zip(sample_tbl["sample"],
                                                            sample_tbl["bam"],
                                                            sample_tbl["condition"])
                       ]

    sys.stderr.write("(source_bam, out_destination) for bams in sample table:\n")


    for src, dest in bam_in_out_paths:
        print(src + " --> " + dest, file=sys.stderr)


    sys.stderr.write("generating softlinks for BAM files...\n")
    [os.symlink(in_bam, out_link_dest) for in_bam, out_link_dest in bam_in_out_paths]

    bai_in_out_paths = [(os.path.abspath(bam + ".bai"),
                        os.path.join(outdir, condition, sample + ".bam.bai")
                        ) for sample, bam, condition in zip(sample_tbl["sample"],
                                                            sample_tbl["bam"],
                                                            sample_tbl["condition"])
                        ]

    sys.stderr.write("(source_bai, out_destination) for bais relative to bam in sample table:\n")

    for src, dest in bai_in_out_paths:
        print(src + " --> " + dest, file=sys.stderr)

    sys.stderr.write("generating softlinks for BAI files...\n")
    [os.symlink(in_bai, out_link_dest) for in_bai, out_link_dest in bai_in_out_paths]

    sys.stderr.write("Done!")
    return None


if __name__ == '__main__':

    description = """Script to soft link all BAMs of the same condition into a specific directory for input to APA-Scan"""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input-sample_table", type=str, dest="sample_tbl",help="path to input sample table generated for APAeval APA-Scan execution workflow run", required=True)
    parser.add_argument("-o", "--outdir", type=str, dest="out_dir", help="path to main output directory under which condition specific subdirectories will be created", required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.sample_tbl, args.out_dir)
