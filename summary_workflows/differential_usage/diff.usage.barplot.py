#!/usr/bin/env python

import os, re, argparse
import pandas as pd
import matplotlib.pyplot as plt


def bar_plot(files, plot_out):
	'''
	Inputs:
		- files: a list of tsv files containing 
			  	gene_id & usage changes (e.g. deltaPPAU in QAPA);
			  	NOTE that these genes have been filtered by FDR/q-value.
			  	File names are like 'AA_diffUTR_03.tsv'
		- [hardcoded] cutoff: the cutoff to determine if the change is big enough or not,
				  e.g. in QAPA, deltaPPAU > 20 means lengthening
		- output: path to store the image
	'''

	files = [f for f in files if '_03.tsv' in f]
	methods = {}
	for file in files:
		method = re.sub('_03.tsv$', '', os.path.basename(file))
		method = re.sub('^.*_', '', method)

		# hardcode the cutoffs since each method has its own calculation
		cutoff = 0
		if method == "QAPA":
			cutoff = 20
		elif method == "diffUTR":
			cutoff = 0.5 #not sure, can ask the author
		# others to be added

		df = pd.read_csv(file, sep='\t', header=None, names=['gene','usage'])
		l = df[df.usage >= cutoff].shape[0]
		s = df[df.usage <= -cutoff].shape[0]
		nc = df.shape[0] - l - s
		methods[method] = [l,s,nc]

	df = pd.DataFrame(methods, index=['Lengthening', 'Shortening', 'NoChange'])
	ax = df.plot.bar(rot=0)
	fig = ax.get_figure()
	fig.savefig(plot_out)


def main(args):
	files = args['tsvfile'] # tsv files
	out_path = args['outdir'][0] # directory to output image
	plot_out = os.path.join(out_path, 'diff.usages.pdf')
	#cutoff = float(args['cutoff'][0])

	bar_plot(files, plot_out)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Differential usage plot')
	#parser.add_argument('--cutoff', required=True, nargs=1)
	parser.add_argument('--tsvfile', required=True, nargs='+')
	parser.add_argument('--outdir', required=True, nargs=1)
	args = vars(parser.parse_args())
	main(args)

# usage: ./diff.usage.barplot.py --tsvfile ../AA_QAPA_03.tsv ../BB_MISO_03.tsv ../AB_diffUTR_03.tsv  --outdir .