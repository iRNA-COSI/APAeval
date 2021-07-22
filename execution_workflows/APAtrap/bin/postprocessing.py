#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Reformat APAtrap identifyDistal3UTR bed file"
    Epilog = "Example usage: python postprocessing.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    return parser.parse_args(args)


def reformat_bed(file_in):
    """
    This function reformats the txt deAPA output file to files for
    identification, quantification, and differential challenges
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    identification_out = open("apatrap_identification_output.bed", "wt")
    quantification_out = open("apatrap_quantification_output.bed", "wt")
    differential_out = open("apatrap_differential_output.tsv", "wt")

    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # replace score column with "."
        # keep just the gene name
        name = row['Gene'].split("|")[0]
        chrom = row['Gene'].split("|")[2]
        score = "."
        strand = row['Gene'].split("|")[0][3]
        apas = row['Predicted_APA'].split(",")

        # write identification file
        for apa in apas:
            chromStart = str(apa)
            # as identified PAS are single-nucleotide, the ending position is the
            # same as starting position
            chromEnd = chromStart
            output = [chrom, chromStart, chromEnd, name, score, strand]
            identification_out.write("\t".join(output) + "\n")

            # write quantification file
            score = str(row['Group_1_1_Total_Exp'])
            output = [chrom, chromStart, chromEnd, name, score, strand]
            quantification_out.write("\t".join(output) + "\n")

        # write differential file
        significance = str(row['p.value'])
        output = [name, significance]
        differential_out.write("\t".join(output) + "\n")

    identification_out.close()
    quantification_out.close()
    differential_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN)


if __name__ == '__main__':
    sys.exit(main())
