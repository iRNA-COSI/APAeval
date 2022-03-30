#!/usr/bin/env python3

import sys
import argparse
import os
from os import listdir
from os.path import isfile, join
import csv

def parse_args(args=None):
    Description = "Reformat CSI-UTR output file to the output file of differential challenge"
    Epilog = "Example usage: python postprocess_differential.py <CSI_UTR_OUTPUT_FOLDER> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("CSI_UTR_OUTPUT_FOLDER", help="CSI-UTR output folder that contains file for differential challenge.")
    parser.add_argument("FILE_OUT", help="Name of output file for the differential challenge.")
    return parser.parse_args(args)


def convert_to_differential(csi_utr_output_folder, file_out):
    """
    This function reformats the txt CSI-UTR output file to files for differential challenge
    :param csi_utr_output_folder: folder containing differential output file to be converted 
    :param file_out: identification challenge output file
    :return: N/A

    example:
    CSI     ENSGENE GENE_SYM        PSI1 (control)  PSI2 (srsf3)    deltaPSI (control-srsf3)        P-value FDR
ENSMUSG00000015652:5736416_5736416-5736328      ENSMUSG00000015652      Steap1  1       1       0       1       1
    """
    files = [f for f in listdir(csi_utr_output_folder) if isfile(join(csi_utr_output_folder, f))]
    # find the file in the format of condition1_VS_condition2_CSI.WITHINUTR.diff.txt"
    for f in files:
       if "diff.txt" in f:
           output_file = join(csi_utr_output_folder, f)
           break

    differential_out = open(file_out, "wt")
    
    with open(output_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        # skip the header line
        next(reader)
        for line in reader:
            gene = line[1]
            pvalue = line[6]
            output = (gene, pvalue)
            # write to output file
            differential_out.write("\t".join(output) + "\n")

    differential_out.close()


def main(args=None):
    args = parse_args(args)
    convert_to_differential(args.CSI_UTR_OUTPUT_FOLDER, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
