#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Reformat DePars output file into the output files of differential challenge"
    Epilog = "Example usage: python convert_output.py <FILE_IN> <FILE_OUT> <MODE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="DaPars output file to be converted.")
    parser.add_argument("FILE_OUT", help="Name of DaPars output for the identification/differential challenge.")
    parser.add_argument("MODE", help="Can either be 'differential' or 'identification")
    return parser.parse_args(args)


def convert_to_differential(file_in, file_out):
    """
    This function reformats the txt deAPA output file to file for
    differential challenges
    :param file_in: txt file to be reformatted
    :param file_out: differential challenge output file
    :return: N/A
    """
    differential_out = open(file_out, "wt")

    df = pd.read_csv(file_in, sep='\t')
    rows = dict()
    for index, row in df.iterrows():
        # write differential file
        # keep just the gene name obtained from the gene column
        # e.g. ENSMUST00000161802.1|NA|chr6|-
        name = row['Gene'].split("|")[0]
        # p val obtained from adjusted.P_val column
        significance = str(row['adjusted.P_val'])
        if name not in rows:
            rows[name] = [significance]
        else:
            rows[name].append(significance)
    # only get the smallest p-value for the corresponding gene
    for gene in rows:
        output = [gene, min(rows[gene])]
        differential_out.write("\t".join(output) + "\n")
    differential_out.close()


def convert_to_identification(file_in, file_out):
    """
    This function reformats the txt deAPA output file to files for
    identification challenge
    :param file_in: txt file to be reformatted
    :param file_out: identification challenge output file
    :return: N/A
    """
    identification_out = open(file_out, "wt")
    identification_outputs = set()
    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # keep just the gene name and orientation from Gene column
        # e.g ENSMUST00000161802.1|NA|chr6|-
        chromosome = row['Gene'].split('|')[2]
        orientation = row['Gene'].split('|')[3]
        proximal_apa = str(row['Predicted_Proximal_APA'])
        # Get distal apa site from the Loci column
        # e.g. chr6:119937615-119938018
        # if orientation is +, distal apa site is the right end of the loci
        if orientation == '+':
            distal_apa = str(row['Loci'].split(':')[1].split('-')[1])
        # else, distal apa site is the left end of the loci
        else:
            distal_apa = str(row['Loci'].split(':')[1].split('-')[0])
        # generate the name for the current identified apa site
        loci = row['Loci'].split(':')[1]
        name = '|'.join([row['Gene'].split('|')[0], chromosome, loci, orientation])

        # add proximal and distal apa sites in separate rows
        identification_outputs.add((chromosome, proximal_apa, proximal_apa, name, '.', orientation))
        identification_outputs.add((chromosome, distal_apa, distal_apa, name, '.', orientation))

    for identification_output in identification_outputs:
        identification_out.write("\t".join(identification_output) + "\n")

    identification_out.close()


def main(args=None):
    args = parse_args(args)
    # check that mode is valid
    if args.MODE == 'identification':
        if not args.FILE_OUT.endswith(".bed"):
            msg = "The identification output file name should end with '.bed'"
            sys.exit(msg)
        convert_to_identification(args.FILE_IN, args.FILE_OUT)
    elif args.MODE == 'differential':
        if not args.FILE_OUT.endswith(".tsv"):
            msg = "The differential output file name should end with '.tsv'"
            sys.exit(msg)
        convert_to_differential(args.FILE_IN, args.FILE_OUT)
    else:
        msg = "The mode set as " + args.MODE +" is not valid. Mode can either be set to 'differential' or 'identification' only"
        sys.exit(msg)


if __name__ == '__main__':
    sys.exit(main())