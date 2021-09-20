#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Create gene symbol file from gft file"
    Epilog = "Example usage: python create_gene_symbol_file.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="GTF file")
    parser.add_argument("FILE_OUT", help="Name of gene symbol text file")
    return parser.parse_args(args)

def create_gene_symbol_file(file_in, file_out):
    """
    This function reads the gtf file and extract information to make
    a gene symbol file
    :param file_in: gtf file to be reformatted
    :return: N/A
    """

    file_out = open(file_out, "wt")
    df = pd.read_csv(file_in, sep='\t', comment="#", header=None)
    for index, row in df.iterrows():
        # read gtf file
        # keep just the gene name and gene symbol
        gene_id = row[8].split(';')[0].split('"')[1]
        gene_name = row[8].split(';')[2].split('"')[1]
        output = [gene_id, gene_name]
        file_out.write("\t".join(output) + "\n")
    file_out.close()

def main(args=None):
    args = parse_args(args)
    create_gene_symbol_file(args.FILE_IN, args.FILE_OUT)

if __name__ == '__main__':
    sys.exit(main())