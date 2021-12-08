#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Reformat APAtrap deAPA output txt file into the output files of differential challenge"
    Epilog = "Example usage: python convert_to_tsv.py <FILE_IN> <DIFFERENTIAL_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("GENOME_FILE", help="GTF genome file for transcript id and gene id information.")
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    parser.add_argument("DIFFERENTIAL_OUT", help="Name of output file for differential challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, differential_out, transcript_id_dict):
    """
    This function reformats the txt deAPA output file to files for
    differential challenges
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    differential_out = open(differential_out, "wt")
    df = pd.read_csv(file_in, sep='\t')
    rows = dict()
    for index, row in df.iterrows():
        # write differential file
        # get the transcript id
        name = row['Gene'].split("|")[0]
        # convert transcript id to gene id
        name = transcript_id_dict[name]
        significance = str(row['p.value'])
        if name not in rows:
            rows[name] = [significance]
        else:
            rows[name].append(significance)
    # only get the smallest p-value for the corresponding gene
    for gene in rows:
        output = [gene, min(rows[gene])]
        differential_out.write("\t".join(output) + "\n")

    differential_out.close()

def create_gene_symbol_file(gtf_file):
    """
    This function reads the gtf file and extract information to make
    a gene symbol file
    :param file_in: gtf file to use to extract transcript id and gene id
    :return: dictionary of transcript id to gene id
    """

    df = pd.read_csv(gtf_file, sep='\t', comment="#", header=None)
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
    return transcript_id_dict


def main(args=None):
    args = parse_args(args)
    transcript_id_dict = create_gene_symbol_file(args.GENOME_FILE)
    reformat_bed(args.FILE_IN, args.DIFFERENTIAL_OUT, transcript_id_dict)


if __name__ == '__main__':
    sys.exit(main())
