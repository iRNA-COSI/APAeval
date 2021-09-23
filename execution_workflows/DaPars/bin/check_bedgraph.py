#!/usr/bin/env python3

import sys
import argparse

def parse_args(args=None):
	Description = "Check DaPars bedgraph input to have leading chr for all sequence regions. Otherwise, throw an error."
	Epilog = "Example usage: python check_bedgraph.py <FILE_IN>"

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument("FILE_IN", help="Input bedgraph file.")
	return parser.parse_args(args)

def check_bedgraph(file_in):
	"""
	This function checks that there is leading chr for all sequence regions
	Otherwise, throw an error
	:param file_in: bedgraph file to be checked
	:return: N/A
	"""
	fin = open(file_in, "rt")

	for line in fin:
		if line[0] != '#' and line[:3] != "chr":
			msg = "Found a row in the sample file " + file_in + " without leading 'chr' in the chromosome column.\n" + \
				" Please use a file with leading 'chr' for all sequence regions."
			sys.exit(msg)
	fin.close()

def main(args=None):
	args = parse_args(args)
	check_bedgraph(args.FILE_IN)


if __name__ == '__main__':
	sys.exit(main())
