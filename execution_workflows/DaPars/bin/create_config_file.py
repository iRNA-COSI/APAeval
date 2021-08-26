#!/usr/bin/env python

import sys
import argparse
import configparser
import os

def parse_args(args=None):
	Description = "Create config file for step 2 of DaPars."
	Epilog = "Example usage: python create_config_file.py" + \
	         " <ANNOTATED_3UTR> <BEDGRAPHS_DIR> <OUTPUT_DIR> <NUM_LEAST_IN_GROUP1> <NUM_LEAST_IN_GROUP2> " \
	         " <COVERAGE_CUTOFF> <FDR_CUTOFF> <PDUI_CUTOFF> <FOLD_CHANGE_CUTOFF> <CONFIG_OUTPUT>"

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument("ANNOTATED_3UTR")
	parser.add_argument("BEDGRAPHS_DIR")
	parser.add_argument("OUTPUT_DIR")
	parser.add_argument("NUM_LEAST_IN_GROUP1")
	parser.add_argument("NUM_LEAST_IN_GROUP2")
	parser.add_argument("COVERAGE_CUTOFF")
	parser.add_argument("FDR_CUTOFF")
	parser.add_argument("PDUI_CUTOFF")
	parser.add_argument("FOLD_CHANGE_CUTOFF")
	parser.add_argument("CONFIG_OUTPUT")

	return parser.parse_args(args)


def get_sample_files(bedgraphs_dir):
	group1 = ""
	group2 = ""
	group = 1
	for folder in os.listdir(bedgraphs_dir):
		if group == 1:
			for file in folder:
				group1 += " " + file
			group += 1
		else:
			for file in folder:
				group2 += " " + file
	return group1, group2


def create_config_file(args):
	"""
	This function creates the config file for step 2 of DaPars
	"""
	config = configparser.ConfigParser()
	config['Annotated_3UTR'] = args.ANNOTATED_3UTR
	group1, group2 = get_sample_files(args.BEDGRAPHS_DIR)
	config['Group1_tophat_aligned_Wig'] = group1
	config['Group2_tophat_aligned_Wig'] = group2
	config['Output_directory'] = args.OUTPUT_DIR
	config['Output_result_file'] = "dapars_output"
	config['Num_least_in_group1'] = args.NUM_LEAST_IN_GROUP1
	config['Num_least_in_group2'] = args.NUM_LEAST_IN_GROUP2
	config['Coverage_cutoff'] = args.COVERAGE_CUTOFF
	config['FDR_cutoff'] = args.FDR_CUTOFF
	config['PDUI_cutoff'] = args.PDUI_CUTOFF
	config['Fold_change_cutoff'] = args.FOLD_CHANGE_CUTOFF
	with open(args.CONFIG_OUTPUT) as configfile:
		config.write(configfile)

def main(args=None):
	args = parse_args(args)
	create_config_file(args)


if __name__ == '__main__':
	sys.exit(main())