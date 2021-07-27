#!/usr/bin/env python3

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
    sample,bamdir,strandinfo,gtf,condition
    """

    input_extensions = []
    sample_info_list = []
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 5
        HEADER = ['sample','bamdir','strandinfo','gtf','condition']
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
            sample, bamdir, strandinfo, gtf, condition = lspl[:len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Sample entry contains spaces!", 'Line', line)
            else:
                print_error("Sample entry has not been specified!", 'Line', line)

            ## Check gtf extension
            if gtf:
                if gtf.find(" ") != -1:
                    print_error("gtf contains spaces!", 'Line', line)
                if not gtf.endswith(".gtf"):
                    print_error("gtf does not have extension '.gtf'", 'Line', line)

            ## Check bamdir entries
            if bamdir:
                if bamdir.find(' ') != -1:
                    print_error("bamdir entry contains spaces!",'Line', line)
                if bamdir.endswith(".bam"):
                    print_error("DO NOT specify the exact BAM file but the directory", 'Line', line)
            
            ## Check strandinfo
            if strandinfo:
                if strandinfo != "reverse_forward" and strandinfo != "unstranded":
                    print_error("Must assign strandinfo to either 'reverse_forward' or 'unstranded'")
            
            ## Create sample mapping dictionary = {group: {replicate : [ sample, bamdir, strandinfo, gtf, condition ]}}
            sample_info = [ sample, bamdir, strandinfo, gtf, condition ]
            sample_info_list.append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_info_list) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(['sample','bamdir','strandinfo','gtf','condition']) + "\n")
            for sample_info in sample_info_list:
                ### Write to file
                fout.write(",".join(sample_info)+"\n")
                

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
