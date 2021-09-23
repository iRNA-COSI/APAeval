#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Reformat DePars output file into the output files of differential challenge"
    Epilog = "Example usage: python convert_output.py <FILE_IN> <MODE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    parser.add_argument("MODE", help="Can either be 'differential' or 'identification")
    return parser.parse_args(args)


def convert_to_differential(file_in):
    """
    This function reformats the txt deAPA output file to file for
    differential challenges
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    differential_out = open("dapars_differential_output.tsv", "wt")

    df = pd.read_csv(file_in, sep='\t')
    rows = dict()
    for index, row in df.iterrows():
        # write differential file
        # keep just the gene name
        name = row['Gene'].split("|")[0]
        significance = str(row['p.value'])
        if name not in rows:
            row[name] = [significance]
        else:
            row[name].append(significance)
    # only get the smallest p-value for the corresponding gene
    for gene in rows:
        output = [gene, min(rows[gene])]
        differential_out.write("\t".join(output) + "\n")
    differential_out.close()


def convert_to_identification(file_in):
    """
    This function reformats the txt deAPA output file to files for
    identification challenge
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    identification_out = open("dapars_identification_output.bed", "wt")
    identification_outputs = set()
    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # keep just the gene name
        chromosome = row['Gene'].split('|')[2]
        proximal_apa = str(row['Predicted_Proximal_APA'])
        orientation = row['Gene'].split('|')[3]
        # if orientation is +, distal apa site is the right end of the loci
        if orientation == '+':
            distal_apa = str(row['Loci'].split(':')[1].split('-')[1])
        # else, distal apa site is the left end of the loci
        else:
            distal_apa = str(row['Loci'].split(':')[1].split('-')[0])
        # the name is the gene id
        loci = row['Loci'].split(':')[1]
        name = '|'.join([row['Gene'].split('|')[0], chromosome, loci, orientation])

        # add proximal and distal apa sites in separate rows
        identification_outputs.add((chromosome, proximal_apa, proximal_apa, name, '.', orientation))
        identification_outputs.add((chromosome, distal_apa, proximal_apa, name, '.', orientation))

    for identification_output in identification_outputs:
        identification_out.write("\t".join(identification_output) + "\n")

    identification_out.close()

def main(args=None):
    args = parse_args(args)
    if args.MODE == 'identification':
        convert_to_identification(args.FILE_IN)
    elif args.MODE == 'differential':
        convert_to_differential(args.FILE_IN)
    else:
        msg = "The mode set as " + args.MODE +" is not valid. Mode can either be set to 'differential' or 'identification' only"
        sys.exit(msg)

if __name__ == '__main__':
    sys.exit(main())