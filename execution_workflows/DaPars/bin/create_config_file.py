#!/usr/bin/env python

import sys
import argparse

def parse_args(args=None):
	Description = "Create config file for step 2 of DaPars."
	Epilog = "Example usage: python create_config_file.py <FILE_IN> <FILE_OUT>"

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

	# The following file is the result of step 1.

	Annotated_3UTR = hg19_refseq_extracted_3UTR.bed

	# A comma-separated list of BedGraph files of samples from condition 1

	Group1_Tophat_aligned_Wig = Condition_A_chrX.wig
	# Group1_Tophat_aligned_Wig=Condition_A_chrX_r1.wig,Condition_A_chrX_r2.wig if multiple files in one group

	# A comma-separated list of BedGraph files of samples from condition 2

	Group2_Tophat_aligned_Wig = Condition_B_chrX.wig

	Output_directory = DaPars_Test_data /

	Output_result_file = DaPars_Test_data

	# At least how many samples passing the coverage threshold in two conditions
	Num_least_in_group1 = 1

	Num_least_in_group2 = 1

	Coverage_cutoff = 30

	# Cutoff for FDR of P-values from Fisher exact test.

	FDR_cutoff = 0.05

	PDUI_cutoff = 0.5

	Fold_change_cutoff = 0.59