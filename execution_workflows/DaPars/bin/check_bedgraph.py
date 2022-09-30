#!/usr/bin/env python3

import sys
import argparse

def parse_args(args=None):
        Description = "Check DaPars bedgraph input to have leading chr for at least one sequence region. Otherwise, throw an error."
        Epilog = "Example usage: python check_bedgraph.py <FILE_IN>"

        parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
        parser.add_argument("FILE_IN", help="Input bedgraph file.")
        return parser.parse_args(args)

def check_bedgraph(file_in):
        """
        This function checks that there is leading chr for at least one sequence region
        Otherwise, throw an error
        :param file_in: bedgraph file to be checked
        :return: N/A
        """
        fin = open(file_in, "rt")
        num_no_chr = 0
        num_line = 0
        for line in fin:
            num_line += 1
            if line[0] != '#' and line[:3] != "chr":
                num_no_chr += 1

        if num_line == num_no_chr:
            msg = "Found all rows in the sample file " + file_in + " without leading 'chr' in the chromosome column. This will cause Dapars to error out.\n" + \
                    " Please use a file with with at least one row with leading 'chr' in the chromosome column."

            sys.exit(msg)

        fin.close()

def main(args=None):
        args = parse_args(args)
        check_bedgraph(args.FILE_IN)


if __name__ == '__main__':
        sys.exit(main())
