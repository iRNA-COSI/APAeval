#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Check DePars input parameters"
    Epilog = "Example usage: python check_input_params.py <FILE_OUT> <MODE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_OUT", help="Name of DaPars output for the identification/differential challenge.")
    parser.add_argument("MODE", help="Can either be 'differential' or 'identification")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)
    # check that mode is valid
    if args.MODE == 'identification':
        if not args.FILE_OUT.endswith(".bed"):
            msg = "The identification output file name should end with '.bed'"
            sys.exit(msg)
    elif args.MODE == 'differential':
        if not args.FILE_OUT.endswith(".tsv"):
            msg = "The differential output file name should end with '.tsv'"
            sys.exit(msg)
    else:
        msg = "The mode set as " + args.MODE +" is not valid. Mode can either be set to 'differential' or 'identification' only"
        sys.exit(msg)


if __name__ == '__main__':
    sys.exit(main())