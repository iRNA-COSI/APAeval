#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Reformat IsoSCM compare step output file into the output file of differential challenge"
    Epilog = "Example usage: python3 postprocess_differential.py <GTF_FILE> <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("GTF_IFLE", help="GTF file to get gene id from.")
    parser.add_argument("FILE_IN", help="IsoSCM compare step output file to be converted.")
    parser.add_argument("FILE_OUT", help="Name of IsoSCM output for differential challenge.")
    return parser.parse_args(args)


def process_gtf(gtf_file):
    gtf = pd.read_csv(gtf_file, sep='\t', comment="#", header=None)
    gene_id_dict = dict()
    interval_list = []
    gene_idx = 0
    for index, row in gtf.itterrows():
        transcript_start = row[3]
        transcript_end = row[4]
        gene_info_list = row[8].split(';')
        gene_id = ""
        for gene_info in gene_info_list:
            if "gene_id" in gene_info.split('"')[0]:
                gene_id = gene_info.split('"')[1]
        if gene_id != "":
            gene_id_dict[gene_idx] = gene_id
            interval_list.append([transcript_start, transcript_end])
            gene_idx += 1

    return gene_id_dict, interval_list


def convert_to_differential(gtf_file, file_in, file_out):
    """
    This function reformats the txt IsoSCM compare ste output file to tsv file for
    differential challenge
    :param file_in: txt file to be reformatted
    :param file_out: differential challenge output file
    :return: N/A
    """
    differential_out = open(file_out, "wt")
    gene_id_dict, interval_list = process_gtf(gtf_file)
    df = pd.read_csv(file_in, sep='\t')
    rows = dict()
    for index, row in df.iterrows():
        # write differential file
        # find gene id of the current change point
        change_point = row['changepoint']
        confidence = row['confidence']
        # find gene id of the current change point
        gene_id = ""
        for i in range(len(interval_list)):
            interval = interval_list[i]
            if interval[0] <= change_point < interval[1]:
                gene_id = gene_id_dict[i]
                break
        if gene_id != "" and gene_id not in rows:
            rows[gene_id] = [confidence]
        else:
            rows[gene_id].append(confidence)
    # only get the smallest p-value for the corresponding gene
    for gene in rows:
        output = [gene, min(rows[gene])]
        differential_out.write("\t".join(output) + "\n")
    differential_out.close()


def main(args=None):
    args = parse_args(args)
    convert_to_differential(args.GTF_FILE, args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())