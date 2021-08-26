#!/usr/bin/env python3

import sys
import argparse

def parse_args(args=None):
	Description = "Reformat DaPars bedgraph input to have leading chr for all sequence regions."
	Epilog = "Example usage: python check_bedgraph.py <FILE_IN> <FILE_OUT>"

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument("FILE_IN", help="Input bedgraph file.")
	parser.add_argument("FILE_OUT", help="Output bedgraph file.")
	return parser.parse_args(args)

def check_bedgraph(file_in, file_out):
	"""
	This function checks that there is leading chr for all sequence regions
	Otherwise, add the leading chr
	:param file_in: bedgraph file to be checked
	:param file_out: output bedgraph file with leading chr
	:return: N/A
	"""
	fin = open(file_in, "rt")
	fout = open(file_out, "wt")

	for line in fin:
		if (line[:3] != "chr"):
			fout.write("chr" + line)
		else:
			fout.write(line)
	fin.close()
	fout.close()

def main(args=None):
	args = parse_args(args)
	check_bedgraph(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
	sys.exit(main())
