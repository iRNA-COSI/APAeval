#!/usr/bin/env python3

import pandas as pd
import sys
import argparse
import os
from os import listdir
from os.path import isfile, join, isdir
import csv


def parse_args(args=None):
    Description = "Reformat CSI-UTR output file to the output file of differential challenge"
    Epilog = "Example usage: python postprocess_differential.py <CSI_UTR_RESULTS_DIR> <RESULT_TYPE> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("CSI_UTR_RESULTS_DIR", help="Path to directory containing CSI-UTR differential expression output results file.")
    parser.add_argument("RESULT_TYPE", choices=["DEXSeq", "PAIRWISE", "WITHIN_UTR"], help="Source/'statistical test type' corresponding to output file. ")
    parser.add_argument("FILE_OUT", help="Name of output file for the differential challenge.")
    return parser.parse_args(args)


def standardize_results_file(df, result_type):
    """
    Process CSI-UTR results tables to a common structure (2 col df of gene_id | p)
    """

    assert result_type in ["DEXSeq", "PAIRWISE", "WITHIN_UTR"]

    # Minimal cols needed for postprocessing - will convert to shared colnames
    dexseq_colmap = {"GENE_ID": "gene_id", "pvalue": "p"}
    pairwise_colmap = {"ENSEMBL_GENE_ID": "gene_id", "P-Value": "p"}
    within_colmap = {"ENSGENE": "gene_id", "P-value": "p"}

    if result_type == "DEXSeq":
        return df.rename(columns=dexseq_colmap)[list(dexseq_colmap.values())]

    elif result_type == "PAIRWISE":
        return df.rename(columns=pairwise_colmap)[list(pairwise_colmap.values())]

    elif result_type == "WITHIN_UTR":
        return df.rename(columns=within_colmap)[list(within_colmap.values())]




def convert_to_differential(csi_utr_results_dir, result_type, file_out):
    """
    This function reformats the txt CSI-UTR output file to files for differential challenge
    :param csi_utr_results_dir: folder containing differential output file to be converted
    :param result_type: string denoting which differential expression approach was used to generate input results files (and will be reported in differential TSV).
    :param file_out: identification challenge output file
    :return: N/A

    example:
    CSI     ENSGENE GENE_SYM        PSI1 (control)  PSI2 (srsf3)    deltaPSI (control-srsf3)        P-value FDR
ENSMUSG00000015652:5736416_5736416-5736328      ENSMUSG00000015652      Steap1  1       1       0       1       1
    """
    # Find results files for specified result type to reprt
    if result_type == "DEXSeq":
        # In DEXSeq results dir, results tables are stored under subdirs (comparison direction)
        subdirs = [join(csi_utr_results_dir, f) for f in listdir(csi_utr_results_dir) if isdir(join(csi_utr_results_dir, f))]

        assert len(subdirs) == 2, f"Expect 2 'DEXSeq' subdirectories within provided directory, {', '.join(subdirs)} were found"

        sys.stderr.write(f"DEXSeq subdirectories found - {', '.join(subdirs)}\n")
        sys.stderr.write(f"P-values are the same since only effect size direction is changed - choosing {subdirs[0]}\n")

        # target file = allCSIs_<comparison>_BYREGION.txt contains results for all input sites
        files = [join(subdirs[0], f) for f in listdir(subdirs[0]) if f.endswith("_BYREGION.txt")]

        # Note that p-values were marginally different with APAeval test date - keeping [0] for consistency with below
        # https://github.com/iRNA-COSI/APAeval/pull/188#issuecomment-1106421434
        results_file = files[0]
        sys.stderr.write(f"Results file found to extract differential challenge outputs - {results_file}\n")

    else:
        # PAIRWISE & WITHIN_UTR have the same structure - just need to find .diff.txt w/in provided dir
        # Again expect two files (with denominator flipped), p-values are the same
        # https://github.com/iRNA-COSI/APAeval/pull/188#issuecomment-1106421434
        files = [join(csi_utr_results_dir, f) for f in listdir(csi_utr_results_dir) if f.endswith(".diff.txt")]

        assert len(files) == 2, f"Expect 2 '*.diff.txt' files in provided directory, following were found - {', '.join(files)}\n"

        sys.stderr.write(f"P-values for two files are the same since only effect size direction is changed -  {files[0]} selected for reporting p-values\n")
        results_file = files[0]

    df = pd.read_csv(results_file, sep='\t')

    # Subset to minimal columns (gene_id & pvalue),
    # renaming to 'gene_id' & 'p' so can uniformly process
    df = standardize_results_file(df, result_type)

    # Pick the smallest p-value per gene_id
    # 'is there any evidence of differential usage within the gene?'
    df = df.groupby("gene_id")["p"].min().reset_index()

    # Output df to headerless TSV of (gene_id | p)
    df.to_csv(file_out, sep="\t", header=False, index=False, na_rep="NA")


def main(args=None):
    args = parse_args(args)
    convert_to_differential(args.CSI_UTR_RESULTS_DIR, args.RESULT_TYPE, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
