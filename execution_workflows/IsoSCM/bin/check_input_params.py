#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Check DePars input parameters"
    Epilog = "Example usage: python check_input_params.py <IDENTIFICATION_OUT_SUFFIX> <RELATIVE_USAGE_QUANTIFICATION_OUT_SUFFIX>" + \
                " <RUN_IDENTIFICATION> <RUN_RELATIVE_USAGE_QUANTIFICATION>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("IDENTIFICATION_OUT_SUFFIX", help="Name of DaPars output for the identification challenge.")
    parser.add_argument("RELATIVE_USAGE_QUANTIFICATION_OUT_SUFFIX", help="Name of DaPars output for the relative usage quantification challenge.")
    parser.add_argument("RUN_IDENTIFICATION", help="Can either be 'true' or 'false'")
    parser.add_argument("RUN_RELATIVE_USAGE_QUANTIFICATION", help="Can either be 'true' or 'false'")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if args.RUN_IDENTIFICATION == True:
        if not args.IDENTIFICATION_OUT_SUFFIX.endswith(".bed"):
            msg = "The identification output file name should end with '.bed'"
            sys.exit(msg) 
    if args.RUN_RELATIVE_USAGE_QUANTIFICATION == True:
        if not args.IDENTIFICATION_OUT_SUFFIX.endswith(".bed"):
            msg = "The relative usage quantification output file name should end with '.bed'"
            sys.exit(msg)


if __name__ == '__main__':
    sys.exit(main())
