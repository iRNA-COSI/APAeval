#!/usr/bin/python3

import sys
import argparse
import csv

def parse_args(args=None):
    Description = "Reformat IsoSCM bed file into the output file of identification challenge"
    Epilog = "Example usage: python postprcoess_identification.py <FILE_IN> <IDENTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input IsoSCM output bed file.")
    parser.add_argument("IDENTIFICATION_OUT", help="Name of output file for identification challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, identification_out):
    """
    This function reformats IsoSCM output file to file for
    identification challenge
    :param file_in: IsoSCM assemble output bed file to be reformatted to identification bed file
    :return: N/A
    """
    identification_out = open(identification_out, "wt")
    with open(file_in) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            chrom = row[0]
            start = row[3]
            end = int(start) + 1
            strand = row[6]
            name = '|'.join([chrom, str(start)+":"+str(end), strand])
            score = "."
            output = (chrom, start, str(end), name, score, strand)
            identification_out.write("\t".join(output) + "\n")

    identification_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.IDENTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
