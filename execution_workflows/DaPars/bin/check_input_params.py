#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Check DePars input parameters"
    Epilog = "Example usage: python check_input_params.py <IDENTIFICATION_OUT_SUFFIX> <DIFFERENTIAL_OUT>" + \
                " <RUN_IDENTIFICATION> <RUN_DIFFERENTIAL>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("IDENTIFICATION_OUT_SUFFIX", help="Name of DaPars output for the identification challenge.")
    parser.add_argument("DIFFERENTIAL_OUT", help="Name of DaPars output for the differential challenge.")
    parser.add_argument("RUN_IDENTIFICATION", help="Can either be 'true' or 'false")
    parser.add_argument("RUN_DIFFERENTIAL", help="Can either be 'true' or 'false")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if args.RUN_IDENTIFICATION == 'true':
        if not args.IDENTIFICATION_OUT_SUFFIX.endswith(".bed"):
            msg = "The identification output file name should end with '.bed'"
            sys.exit(msg)
    elif args.RUN_DIFFERENTIAL == 'true':
        if not args.DIFFERENTIAL_OUT.endswith(".tsv"):
            msg = "The differential output file name should end with '.tsv'"
            sys.exit(msg)


if __name__ == '__main__':
    sys.exit(main())