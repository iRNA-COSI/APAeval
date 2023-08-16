#!/usr/bin/env python3

import sys
import argparse
import configparser
import os

def parse_args(args=None):
	Description = "Create config file for step 2 of DaPars."
	Epilog = "Example usage: python create_config_file.py" + \
			" <ANNOTATED_3UTR> <BEDGRAPHS_DIR> <OUTPUT_DIR> <NUM_LEAST_IN_GROUP1> <NUM_LEAST_IN_GROUP2>"  + \
			" <COVERAGE_CUTOFF> <FDR_CUTOFF> <PDUI_CUTOFF> <FOLD_CHANGE_CUTOFF> <CONFIG_OUTPUT>" + \
			" <BEDGRAPH_FILE> <RUN_MODE>"

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
	parser.add_argument("BEDGRAPH_FILE")
	parser.add_argument("RUN_MODE")

	return parser.parse_args(args)


def get_sample_files_differential(bedgraphs_dir):
	group1 = []
	group2 = []
	group = 1
	for folder in os.listdir(bedgraphs_dir):
		if group == 1:
			for file in os.listdir(os.path.join(bedgraphs_dir,folder)):
				group1.append(os.path.join(bedgraphs_dir, folder, file))
			group += 1
		else:
			for file in os.listdir(os.path.join(bedgraphs_dir,folder)):
				group2.append(os.path.join(bedgraphs_dir, folder, file))
	group1 = ','.join(group1)
	group2 = ','.join(group2)
	return group1, group2


def get_sample_files(bedgraph_dir, bedgraph_file, run_mode):
	if run_mode == "identification" or run_mode == 'relative_usage_quantification':
		return bedgraph_file, bedgraph_file
	else:
		return get_sample_files_differential(bedgraph_dir)


def create_config_file(args):
	"""
	This function creates the config file for step 2 of DaPars
	"""
	group1, group2 = get_sample_files(args.BEDGRAPHS_DIR, args.BEDGRAPH_FILE, args.RUN_MODE)
	config = {
		'Annotated_3UTR': args.ANNOTATED_3UTR,
		'Group1_Tophat_aligned_Wig': group1,
		'Group2_Tophat_aligned_Wig': group2,
		'Output_result_file': 'dapars_output',
		'Output_directory': args.OUTPUT_DIR,
		'Num_least_in_group1': args.NUM_LEAST_IN_GROUP1,
		'Num_least_in_group2': args.NUM_LEAST_IN_GROUP2,
		'Coverage_cutoff': args.COVERAGE_CUTOFF,
		'FDR_cutoff': args.FDR_CUTOFF,
		'PDUI_cutoff': args.PDUI_CUTOFF,
		'Fold_change_cutoff': args.FOLD_CHANGE_CUTOFF
	}

	with open(args.CONFIG_OUTPUT, 'w') as f:
		for key, value in config.items():
			f.write('%s=%s\n' % (key, value))


def main(args=None):
	args = parse_args(args)
	create_config_file(args)


if __name__ == '__main__':
	sys.exit(main())
