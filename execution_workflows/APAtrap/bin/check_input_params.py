#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

"""
This file checks if mode and the extension of the file names set in conf/modules.config are valid
"""


def parse_args(args=None):
    Description = "Check DePars input parameters"
    Epilog = "Example usage: python check_input_params.py <IDENTIFICATION_OUT_SUFFIX> <QUANTIFICATION_OUT_SUFFIX>" +\
             " <DIFFERENTIAL_OUT> <RUN_IDENTIFICATION> <RUN_QUANTIFICATION> <RUN_DIFFERENTIAL>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("IDENTIFICATION_OUT_SUFFIX", help="Name of APAtrap output for the identification challenge.")
    parser.add_argument("QUANTIFICATION_OUT_SUFFIX", help="Name of APAtrap output for the quantification challenge.")
    parser.add_argument("DIFFERENTIAL_OUT", help="Name of APAtrap output for the differential challenge.")
    parser.add_argument("RUN_IDENTIFICATION", help="Can either be 'true' or 'false")
    parser.add_argument("RUN_QUANTIFICATION", help="Can either be 'true' or 'false")
    parser.add_argument("RUN_DIFFERENTIAL", help="Can either be 'true' or 'false")
    return parser.parse_args(args)


def main(args=None):
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
    elif not args.RUN_DIFFERENTIAL:
        if not args.FILE_OUT.endswith(".tsv"):
            msg = "The differential output file name should end with '.tsv'"
            sys.exit(msg)
        # check that there are


if __name__ == '__main__':
    sys.exit(main())