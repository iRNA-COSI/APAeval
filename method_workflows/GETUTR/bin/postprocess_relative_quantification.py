#!/usr/bin/python2

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    description = "Reformat GETUTR output bed file into the output file of relative usage quantification challenges"
    epilog = "Example usage: python convert_to_bed.py <FILE_IN> <RELATIVE_USAGE_QUANTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("FILE_IN", help="Input GETUTR output txt file.")
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
    
    columns = ['chr', 'start', 'end', 'transcript', 'score', 'strand']
    df = pd.read_csv(file_in, sep='\t', header=None, names=columns)

    for transcript in df['transcript'].unique():
        count = 0
        for index, row in df[df.transcript == transcript].sort_values(by='start').iterrows():
            relative_usage_quantification_outputs.append([row.chr, row.start, row.end, transcript + "_" + str(count), row.score, row.strand])
            count += 1
    ret_df = pd.DataFrame(relative_usage_quantification_outputs)
    ret_df.to_csv(relative_usage_quantification_out, sep='\t', header=False, index=False)


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.RELATIVE_USAGE_QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
