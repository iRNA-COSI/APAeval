#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Reformat DePars output file into the output files of differential challenge"
    Epilog = "Example usage: python convert_to_tsv.py <FILE_IN> <IDENTIFICATION_OUT> <DIFFERENTIAL_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    parser.add_argument("IDENTIFICATION_OUT", help="Name of output file for identification challenge")
    parser.add_argument("DIFFERENTIAL_OUT", help="Name of output file for differential challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, differential_out):
    """
    This function reformats the txt deAPA output file to files for
    differential challenges
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    differential_out = open(differential_out, "wt")

    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # write differential file
        # keep just the gene name
        name = row['Gene'].split("|")[0]
        significance = str(row['p.value'])
        output = [name, significance]
        differential_out.write("\t".join(output) + "\n")

    differential_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())