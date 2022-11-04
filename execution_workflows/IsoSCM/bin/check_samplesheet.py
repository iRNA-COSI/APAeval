#!/usr/bin/env python

import os
import sys
import argparse

def parse_args(args=None):
    Description = "Reformat APAeval samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context='Line', context_str=''):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != '' and context_str != '':
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(error, context.strip(), context_str.strip())
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample, bam, strand
    """

    input_extensions = []
    sample_info_list = []
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 4
        HEADER = ['sample', 'bam', 'strand', 'read_type']
        header = fin.readline().strip().split(",")
        if header[:len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip() for x in line.strip().split(",")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error("Invalid number of columns (minimum = {})!".format(len(HEADER)), 'Line', line)

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error("Invalid number of populated columns (minimum = {})!".format(MIN_COLS), 'Line', line)

            ## Check group name entries
            sample, bam, strand, read_type = lspl[:len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Sample entry contains spaces!", 'Line', line)
            else:
                print_error("Sample entry has not been specified!", 'Line', line)

            ## Check bam extension
            if bam:
                if bam.find(" ") != -1:
                    print_error("bam contains spaces!", 'Line', line)
                if not bam.endswith(".bam"):
                    print_error("bam does not have extension 'bam'", 'Line', line)
            
            ## Check the strand
            if strand:
                if strand != "reverse_forward" and strand != "unstranded":
                    print_error("strand can only be either reverse_forward or unstranded, received " + strand, 'Line', line)
            
            ## Check the read_type
            if read_type:
                if read_type != "single" and read_type != "paired":
                    print_error("read_type can only be either single or paired, received " + read_type, 'Line', line)

            # Check strand and read_type combination
            if strand == "reverse_forward" and read_type == "single":
                print_error("For single end read, strand should be set to unstranded to avoid error thrown by the tool")

            ## Create sample mapping dictionary = {group: {replicate : [ barcode, input_file, genome, gtf, is_transcripts ]}}
            sample_info = [ sample, bam, strand, read_type]
            sample_info_list.append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_info_list) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(['sample','bam', 'strand', 'read_type']) + "\n")
            for sample_info in sample_info_list:
                ### Write to file
                fout.write(",".join(sample_info)+"\n")

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())

