#!/usr/bin/env python3

import sys
import argparse
import os.path
import pandas as pd


def parse_args(args=None):
    Description = "Create sample information file"
    Epilog = "Example usage: python create_sample_information_file.py <SAMPLESHEET> <SAMPLE_INFORMATION_FILE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("SAMPLESHEET", help="Path to samplesheet")
    parser.add_argument("SAMPLE_INFORMATION_FILE", help="Path to sample information file")
    return parser.parse_args(args)


def write_to_sample_information_file(args):
    """
        This function writes to sampleInformation.txt
        :param args: containing samplesheet for sampleInformation.txt
        :return: N/A
    """
    outputs = []
    df = pd.read_csv(args.SAMPLESHEET)
    # default value for number of mapped reads as specified by CSI-UTR's readme
    num_mapped_reads = 1000000
    for index, row in df.iterrows():
        # bam file name of the sample
        bam_file_name = os.path.basename(row['bam'])
        # condition name of the sample
        condition = row['condition']
        # replicate number of the sample
        replicate = row['replicate']
        output = [bam_file_name, condition, replicate, num_mapped_reads]
        outputs.append(output)
    df = pd.DataFrame(outputs)
    df.to_csv(args.SAMPLE_INFORMATION_FILE, sep='\t', header=False, index=False)


def main(args=None):
    args = parse_args(args)
    write_to_sample_information_file(args)


if __name__ == '__main__':
    sys.exit(main())