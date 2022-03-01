#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    description = "Reformat APAtrap identifyDistal3UTR bed file into the output files of " + \
                  "identification and quantification challenges"
    epilog = "Example usage: python convert_to_tsv.py <FILE_IN> <RUN_IDENTIFICATION> " + \
             "<RUN_QUANTIFICATION> <IDENTIFICATION_OUT> <QUANTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    parser.add_argument("RUN_IDENTIFICATION", help="Boolean indicating whether identification run mode is on")
    parser.add_argument("RUN_QUANTIFICATION", help="Boolean indicating whether quantification run mode is on")
    parser.add_argument("IDENTIFICATION_OUT", help="Name of output file for identification challenge")
    parser.add_argument("QUANTIFICATION_OUT", help="Name of output file for quantification challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, run_identification, run_quantification, identification_out, quantification_out):
    """
    This function reformats the txt predictAPA output file to files for
    identification and quantification challenges
    :param file_in: txt file to be reformatted
    :param run_identification: boolean whether to write identification file
    :param run_quantification: boolean whether to write quantification file
    :param identification_out: identification output file name
    :param quantification_out: quantification output file name
    :return: N/A
    """

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
        proximal_apa_sites = row['Predicted_APA'].split(",")

        # Write identification file
        for proximal_apa_site in proximal_apa_sites:
            proximal_apa_site = int(proximal_apa_site)
            # replace score column with "."
            score = "."
            # As identified PAS are single-nucleotide, the ending position is the
            # starting position + 1 since bed file is end exclusive
            output = (chrom, str(proximal_apa_site), str(proximal_apa_site + 1), name, score, strand)
            identification_outputs.add(output)

        # Get distal apa site from the Loci column
        # e.g. chr6:119937615-119938018
        # if orientation is +, distal apa site is the right end of the loci
        if strand == '+':
            # -1 because bed file coordinate is 0 based and end exclusive
            distal_apa_site = int(row['Loci'].split(':')[1].split('-')[1]) - 1
            # else, distal apa site is the left end of the loci
        else:
            distal_apa_site = int(row['Loci'].split(':')[1].split('-')[0])
            # generate the name for the current identified apa site
        output = (chrom, str(distal_apa_site), str(distal_apa_site + 1), name, score, strand)
        identification_outputs.add(output)

        # Write quantification file
        all_apa_sites = [int(x) for x in proximal_apa_sites]
        all_apa_sites.append(distal_apa_site)

        # Expression level of each APA sites (from the most proximal site to the distal site, seperated by comma).
        expressions = row['Group_1_1_Separate_Exp'].split(',')

        for i in range(len(all_apa_sites)):
            output = (chrom, str(all_apa_sites[i]), str(all_apa_sites[i] + 1), name, str(expressions[i]), strand)
            quantification_outputs.add(output)

    # write to files
    if run_identification == 'true':
        identification_out = open(identification_out, "wt")
        for row in identification_outputs:
            identification_out.write("\t".join(row) + "\n")
        identification_out.close()
    if run_quantification == 'true':
        quantification_out = open(quantification_out, "wt")
        for row in quantification_outputs:
            quantification_out.write("\t".join(row) + "\n")
        quantification_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.RUN_IDENTIFICATION, args.RUN_QUANTIFICATION,
                 args.IDENTIFICATION_OUT, args.QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
