import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bed', help='The BED file containing predictions merged with ground truth.')
parser.add_argument('-o', help='Name of the output json file.')
args = parser.parse_args()

fname = args.bed
# output file name
if not args.o:
    oname = fname[:-4] + '_corr.json'
else:
    oname = args.o

out = pd.read_csv(fname, delimiter='\t')

# initialize vectors for ground truth sites and prediction
vec_true = []
vec_pred = []

# multiple predicted sites for one ground truth?
multiple_predicted_sites = out.duplicated(['chrom_g', 'chromStart_g', 'chromEnd_g', 'strand_g'], keep=False)

# iterate over matched sites
for (idx, row), is_mult in zip(out.iterrows(), multiple_predicted_sites):
    if is_mult:
        # sum up all prediction sites that were assigned to this ground truth site
        # needs to be implemented or can be skipped for now since these are usually only a few
        pass
    elif row['score_g'] == 0: # if there was no ground truth match, expression was set to 0 and the site is excluded
        pass
    else:
        vec_true.append(row['score_g'])
        # weighted expression in case there are multiple ground truth sites for one predicted site
        vec_pred.append(row['score_p']*row['weight'])

# correlation coefficient
r = pearsonr(vec_true, vec_pred)[0]

# export json
with open(oname, 'w') as f:
    json.dump({'correlation_coefficient': r}, f)

print('correlation coefficient: %.2f'%(r))
