#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Reformat APAtrap identifyDistal3UTR bed file into the output files of identification and quantification challenges"
    Epilog = "Example usage: python convert_to_tsv.py <FILE_IN> <IDENTIFICATION_OUT> <QUANTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    parser.add_argument("IDENTIFICATION_OUT", help="Name of output file for identification challenge")
    parser.add_argument("QUANTIFICATION_OUT", help="Name of output file for quantification challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, identification_out, quantification_out):
    """
    This function reformats the txt predictAPA output file to files for
    identification and quantification challenges
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    identification_out = open(identification_out, "wt")
    quantification_out = open(quantification_out, "wt")

    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # replace score column with "."
        # keep just the gene name
        name = row['Gene'].split("|")[0]
        chrom = row['Gene'].split("|")[2]
        strand = row['Gene'].split("|")[3]
        apas = row['Predicted_APA'].split(",")

        # write identification file
        for apa in apas:
            chromStart = str(apa)
            # as identified PAS are single-nucleotide, the ending position is the
            # same as starting position
            chromEnd = chromStart
            score = "."
            output = [chrom, chromStart, chromEnd, name, score, strand]
            identification_out.write("\t".join(output) + "\n")

            # write quantification file
            score = str(row['Group_1_1_Total_Exp'])
            output = [chrom, chromStart, chromEnd, name, score, strand]
            quantification_out.write("\t".join(output) + "\n")

    identification_out.close()
    quantification_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.IDENTIFICATION_OUT, args.QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
