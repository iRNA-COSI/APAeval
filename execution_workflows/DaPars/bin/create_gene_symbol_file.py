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
    transcript_id_dict = dict()
    for index, row in df.iterrows():
        # read gtf file
        # keep just the gene id and gene name
        gene_info_list = row[8].split(';')
        transcript_id = ""
        for gene_info in gene_info_list:
            if "transcript_id" in gene_info.split('"')[0]:
                transcript_id = gene_info.split('"')[1]
            if "gene_id" in gene_info.split('"')[0]:
                gene_id = gene_info.split('"')[1]
        if transcript_id != "":
            transcript_id_dict[transcript_id] = gene_id
    file_out.write("# transcript id      gene id" + "\n")
    for transcript_id in transcript_id_dict:
        output = [transcript_id, transcript_id_dict[transcript_id]]
        file_out.write("\t".join(output) + "\n")
    file_out.close()

def main(args=None):
    args = parse_args(args)
    create_gene_symbol_file(args.FILE_IN, args.FILE_OUT)

if __name__ == '__main__':
    sys.exit(main())