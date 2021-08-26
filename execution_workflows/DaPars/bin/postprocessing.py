#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Reformat DaPars output file"
    Epilog = "Example usage: python postprocessing.py <FILE_IN>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input DaPars output file.")
    return parser.parse_args(args)


def reformat_bed(file_in):
    """
    This function reformats the DaPars output file to files for
    differential challenge
    :param file_in: file to be reformatted
    :return: N/A
    """
    differential_out = open("dapars_differential_output.tsv", "wt")

    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # replace score column with "."
        # keep just the gene name
        name = row['Gene'].split("|")[0]

        # write differential file
        significance = str(row['P_val'])
        output = [name, significance]
        differential_out.write("\t".join(output) + "\n")

    differential_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN)


if __name__ == '__main__':
    sys.exit(main())
