#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

"""
This file checks if the extensions of the file names set in conf/modules.config are valid for each mode
and checks that if run_differential is true, there should be exactly two distinct conditions in 
the input sample file
"""


def parse_args(args=None):
    Description = "Check DePars input parameters"
    Epilog = "Example usage: python check_input_params.py <SAMPLE_FILE> <IDENTIFICATION_OUT_SUFFIX> <QUANTIFICATION_OUT_SUFFIX>" +\
             " <DIFFERENTIAL_OUT> <RUN_IDENTIFICATION> <RUN_QUANTIFICATION> <RUN_DIFFERENTIAL>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("SAMPLE_FILE", help="Input sample file for the run.")
    parser.add_argument("IDENTIFICATION_OUT_SUFFIX", help="Name of APAtrap output for the identification challenge.")
    parser.add_argument("QUANTIFICATION_OUT_SUFFIX", help="Name of APAtrap output for the quantification challenge.")
    parser.add_argument("RELATIVE_USAGE_QUANTIFICATION_OUT_SUFFIX", help="Name of APAtrap output for the relative usage quantification challenge.")
    parser.add_argument("DIFFERENTIAL_OUT", help="Name of APAtrap output for the differential challenge.")
    parser.add_argument("RUN_IDENTIFICATION", help="Can either be 'true' or 'false")
    parser.add_argument("RUN_QUANTIFICATION", help="Can either be 'true' or 'false")
    parser.add_argument("RUN_RELATIVE_USAGE_QUANTIFICATION", help="Can either be 'true' or 'false")
    parser.add_argument("RUN_DIFFERENTIAL", help="Can either be 'true' or 'false")
    return parser.parse_args(args)


def main(args=None):
    # parse arguments
    args = parse_args(args)

    # check that output file names are valid
    if args.RUN_IDENTIFICATION:
        if not args.IDENTIFICATION_OUT_SUFFIX.endswith(".bed"):
            msg = "The identification output file name should end with '.bed'"
            sys.exit(msg)
    if args.RUN_QUANTIFICATION:
        if not args.QUANTIFICATION_OUT_SUFFIX.endswith(".bed"):
            msg = "The quantification output file name should end with '.bed'"
            sys.exit(msg)
    if args.RUN_RELATIVE_USAGE_QUANTIFICATION:
        if not args.RELATIVE_USAGE_QUANTIFICATION_OUT_SUFFIX.endswith(".bed"):
            msg = "The relative usage quantification output file name should end with '.bed'"
            sys.exit(msg)
    if args.RUN_DIFFERENTIAL:
        if not args.DIFFERENTIAL_OUT.endswith(".tsv"):
            msg = "The differential output file name should end with '.tsv'"
            sys.exit(msg)
        # check that there are exactly two distinct conditions
        sample = pd.read_csv(args.SAMPLE_FILE)
        num_conditions = len(set(sample["condition"]))
        if num_conditions != 2:
            msg = "The number of distinct conditions in the input sample file to run differential should be 2.\n"
            msg += "Got " + str(num_conditions) + "."
            sys.exit(msg)


if __name__ == '__main__':
    sys.exit(main())
