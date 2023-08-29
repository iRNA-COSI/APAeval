#!/usr/bin/env python3

import sys
import argparse
import os
import csv

def parse_args(args=None):
    Description = "Reformat CSI-UTR output file to the output file of identification challenge"
    Epilog = "Example usage: python postprocess_identification.py <FILE_IN> <BAM> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="CSI-UTR output file to be converted.")
    parser.add_argument("BAM", help="Name of input BAM file to CSI-UTR")
    parser.add_argument("FILE_OUT", help="Suffix of CSI-UTR output for the identification challenge.")
    return parser.parse_args(args)


def convert_to_identification(file_in, bam, file_out):
    """
    This function reformats the txt CSI-UTR output file to files for identification challenge
    :param file_in: txt file to be reformatted
    :param file_out: identification challenge output file
    :return: N/A

    example:
    chr1    134200409       134201252       ENSMUSG00000042451:134200409_134200409-134201252        1       +       0
    chr1    134201356       134201604       ENSMUSG00000042429:134202951_134201604-134201356        1       -       0
    """
    bam_name = os.path.basename(bam)
    bam_name = os.path.splitext(bam_name)[0]
    file_in += bam_name + ".CSIcoverage"

    identification_out = open(file_out, "wt")
    
    with open(file_in, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            chromosome = line[0]
            strand = line[5]
            # the third column contains the apa site
            # if orientation is positive, this is the end of the CSI site
            # do -1 for the start of apa site since bed files are end exclusive 
            if strand == '+':
                end = line[2]
                start = int(end) - 1
            # if orientation is negative, the apa site is the start of the CSI site
            else:
                start = line[2]
                end = int(start) + 1
            score = "."
            loci = str(start) + "-" + str(end)
            name = '|'.join([chromosome, loci, strand])
            output = (chromosome, str(start), str(end), name, score, strand)
            
            # write to output file
            identification_out.write("\t".join(output) + "\n")

    identification_out.close()


def main(args=None):
    args = parse_args(args)
    convert_to_identification(args.FILE_IN, args.BAM, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
