#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Reformat DePars output file into the output files of differential challenge"
    Epilog = "Example usage: python convert_output.py <FILE_IN> <FILE_OUT> <MODE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="DaPars output file to be converted.")
    parser.add_argument("FILE_OUT", help="Name of DaPars output for the identification/relative usage quantification/differential challenge.")
    parser.add_argument("MODE", help="Can either be 'identification', 'relative_usage_quantification', or 'differential'")
    return parser.parse_args(args)


def parse_deAPA_output(file_in):
    outputs = set()
    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # keep just the gene name and orientation from Gene column
        # e.g ENSMUST00000203335.1|ENSMUSG00000045962.16|chr6|-
        chromosome = row['Gene'].split('|')[2]
        orientation = row['Gene'].split('|')[3]
        proximal_apa = int(row['Predicted_Proximal_APA'])
        distal_relative_usage = float(row['A_1_PDUI'])
        proximal_relative_usage = 1 - distal_relative_usage

        # Get distal apa site from the Loci column
        # e.g. chr6:119937615-119938018
        # if orientation is +, distal apa site is the right end of the loci
        if orientation == '+':
            # -1 because bed file coordinate is 0 based and end exclusive
            distal_apa = int(row['Loci'].split(':')[1].split('-')[1]) - 1
        # else, distal apa site is the left end of the loci
        else:
            #TODO: remove -1 from the line below once Dapars has fixed the error of reporting 1 nt downstream of the actual annotated start
            distal_apa = int(row['Loci'].split(':')[1].split('-')[0]) - 1
        # generate the name for the current identified apa site
        loci = row['Loci'].split(':')[1]
        name = '|'.join([row['Gene'].split('|')[0], chromosome, loci, orientation])

        # add proximal and distal apa sites in separate rows
        outputs.add((chromosome, str(proximal_apa), str(proximal_apa+1), name, str(proximal_relative_usage), orientation))
        outputs.add((chromosome, str(distal_apa), str(distal_apa+1), name, str(distal_relative_usage), orientation))

    return outputs
        

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
        # keep just the gene id obtained from the gene column
        # this column has transcript id, gene id, chromosome, orientation
        # e.g. ENSMUST00000203335.1|ENSMUSG00000045962.16|chr6|-
        name = row['Gene'].split("|")[1]
        # p val obtained from P_val column
        significance = row['P_val']
        if name not in rows:
            rows[name] = [significance]
        else:
            rows[name].append(significance)
    # only get the smallest p-value for the corresponding gene
    for gene in rows:
        output = [gene, str(min(rows[gene]))]
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
    identification_outputs = parse_deAPA_output(file_in)
    for identification_output in identification_outputs:
        # change the score column to '.'
        identification_output = list(identification_output)
        identification_output[4] = "."
        identification_out.write("\t".join(identification_output) + "\n")

    identification_out.close()


def convert_to_relative_usage_quantification(file_in, file_out):
    """
    This function reformats the txt deAPA output file to files for
    relative usage quantification challenge
    :param file_in: txt file to be reformatted
    :param file_out: relative usage quantification challenge output file
    :return: N/A
    """
    relative_usage_quantification_out = open(file_out, "wt")
    relative_usage_quantification_outputs = parse_deAPA_output(file_in)
    for relative_usage_quantification_output in relative_usage_quantification_outputs:
        relative_usage_quantification_out.write("\t".join(relative_usage_quantification_output) + "\n")

    relative_usage_quantification_out.close()


def main(args=None):
    args = parse_args(args)
    # check that mode is valid
    if args.MODE == 'identification':
        convert_to_identification(args.FILE_IN, args.FILE_OUT)
    elif args.MODE == 'relative_usage_quantification':
       convert_to_relative_usage_quantification(args.FILE_IN, args.FILE_OUT)
    elif args.MODE == 'differential':
        convert_to_differential(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
