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

    identification_outputs = set()
    quantification_outputs = set()

    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        # Get name, chromosome, strand, from the Gene column
        # e.g.ENSMUST00000023572.14|NA|16|+
        # Generate name for the apa
        id = row['Gene'].split("|")[0]
        chrom = row['Gene'].split("|")[2]
        strand = row['Gene'].split("|")[3]
        loci = row['Loci'].split(':')[1]
        name = '|'.join([id, chrom, loci, strand])

        # column with proximal apa sites
        # in decreasing order when strand is - e.g 119924739,119924260
        # in increasing order when strand i + e.g. 78340249,78340414
        apas = row['Predicted_APA'].split(",")

        # Write identification file
        for apa in apas:
            proximal_apa = int(apa)
            # replace score column with "."
            score = "."
            # As identified PAS are single-nucleotide, the ending position is the
            # same as starting position, end exclusive
            output = (chrom, str(proximal_apa), str(proximal_apa + 1), name, score, strand)
            identification_outputs.add(output)

        # Get distal apa site from the Loci column
        # e.g. chr6:119937615-119938018
        # if orientation is +, distal apa site is the right end of the loci
        if strand == '+':
            # -1 because bed file coordinate is 0 based and end exclusive
            distal_apa = int(row['Loci'].split(':')[1].split('-')[1]) - 1
            # else, distal apa site is the left end of the loci
        else:
            distal_apa = int(row['Loci'].split(':')[1].split('-')[0])
            # generate the name for the current identified apa site
        output = (chrom, str(distal_apa), str(distal_apa + 1), name, score, strand)
        identification_outputs.add(output)

        # Write quantification file
        apas = [int(x) for x in apas]
        apas.append(distal_apa)

        # Expression level of each APA sites (from the most proximal site to the distal site, seperated by comma).
        expressions = row['Group_1_1_Separate_Exp'].split(',')

        for i in range(len(apas)):
            output = (chrom, str(apas[i]), str(apas[i] + 1), name, str(expressions[i]), strand)
            quantification_outputs.add(output)

    # write to files
    for row in identification_outputs:
        identification_out.write("\t".join(row) + "\n")
    for row in quantification_outputs:
        quantification_out.write("\t".join(row) + "\n")

    identification_out.close()
    quantification_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.IDENTIFICATION_OUT, args.QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
