#!/usr/bin/env python

import os, re, argparse
import pandas as pd
import matplotlib.pyplot as plt


def bar_plot(df_path, cutoff, plot_out):
	'''
	Inputs:
		- df_path: path to a directory of tsv files containing 
			  	gene_id & usage changes (e.g. deltaPPAU in QAPA);
			  	NOTE that these genes have been filtered by FDR/q-value.
			  	File names are like 'AA_diffUTR_03.tsv'
		- cutoff: the cutoff to determine if the change is big enough or not,
				  e.g. in QAPA, deltaPPAU > 20 means lengthening
		- output: path to store the image
	'''

	files = os.listdir(df_path)
	files = [f for f in files if '_03.tsv' in f]

	methods = {}
	for file in files:
		method = re.sub('_03.tsv$', '', file)
		method = re.sub('^.*_', '', method)
		df = pd.read_csv(os.path.join(df_path, file), sep='\t', header=None, names=['gene','usage'])
		l = df[df.usage >= cutoff].shape[0]
		s = df[df.usage <= -cutoff].shape[0]
		nc = df.shape[0] - l - s
		methods[method] = [l,s,nc]

	df = pd.DataFrame(methods, index=['Lengthening', 'Shortening', 'NoChange'])
	ax = df.plot.bar(rot=0)
	fig = ax.get_figure()
	fig.savefig(plot_out)


def main(args):
	df_path = args['pathToTSV'][0] # directory to tsv files
	out_path = args['out'][0] # directory to output image
	plot_out = out_path + '/diff.usages.pdf'
	cutoff = float(args['cutoff'][0])

	bar_plot(df_path,cutoff,plot_out)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Differential usage plot')
	parser.add_argument('--cutoff', required=True, nargs=1)
	parser.add_argument('--pathToTSV', required=True, nargs=1)
	parser.add_argument('--out', required=True, nargs=1)
	args = vars(parser.parse_args())
	main(args)

# usage: ./diff.usage.barplot.py --pathToTSV "../" --out "." --cutoff 0.5