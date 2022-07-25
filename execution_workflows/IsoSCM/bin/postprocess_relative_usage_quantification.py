#!/usr/bin/python3

import sys
import argparse
import csv

def parse_args(args=None):
    Description = "Reformat IsoSCM bed file into the output file of identification challenge"
    Epilog = "Example usage: python postprcoess_identification.py <FILE_IN> <RUN_IDENTIFICATION> <RELATIVE_USAGE_QUANTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input IsoSCM output bed file.")
    parser.add_argument("RELATIVE_USAGE_QUANTIFICATION_OUT", help="Name of otuput file for relative usage quantification challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, relative_usage_quantification_out):
    """
    This function reformats IsoSCM output file to file for
    relative usage quantification challenge
    :param file_in: IsoSCM compare step output file to be reformatted to relative usage quantification bed file
    :return: N/A
    """
    relative_usage_quantification_out = open(relative_usage_quantification_out, "wt")
    with open(file_in) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            chrom = row[0]
            # need to do -1 to convert from 1-based gtf file coordinate to 0-based bed file coordinate
            start = int(row[3]) - 1
            end = row[3]
            strand = row[6]
            name = '|'.join([chrom, str(start) + ":" + str(end), strand])
            score = "."
            output = (chrom, str(start), str(end), name, score, strand)
            relative_usage_quantification_out.write("\t".join(output) + "\n")

    relative_usage_quantification_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.RELATIVE_USAGE_QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
