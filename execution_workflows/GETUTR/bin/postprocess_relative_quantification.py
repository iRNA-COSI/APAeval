#!/usr/bin/env python

import sys
import argparse


def parse_args(args=None):
    description = "Reformat GETUTR output bed file into the output file of relative usage quantification challenges"
    epilog = "Example usage: python convert_to_bed.py <FILE_IN> <RUN_IDENTIFICATION> <RELATIVE_USAGE_QUANTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("FILE_IN", help="Input deAPA output txt file.")
    parser.add_argument("RELATIVE_USAGE_QUANTIFICATION_OUT", help="Name of output file for relative usage quantification challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, relative_usage_quantification_out):
    """
    This function reformats GETUTR bed output file to file for relative quantification challenges
    :param file_in: bed file to be reformatted
    :param relative_usage_quantification_out: relative usage quantification output file name
    :return: N/A
    """

    relative_usage_quantification_outputs = []
    
    outfile = open(relative_usage_quantification_out, "w")

    with open(file_in) as f:
        curr_transcript = ''
        count = 0
        for line in f:
            chrom = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            transcript = line.split('\t')[3]
            if curr_transcript == transcript:
                count += 1
            else:
                curr_transcript = transcript
                count = 0
            name = transcript + '_' + str(count)
            score = line.split('\t')[4]
            strand = line.split('\t')[5]
            row = [chrom, start, end, name, score, strand]
            outfile.write("\t".join(row))
    outfile.close()

def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.RELATIVE_USAGE_QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
