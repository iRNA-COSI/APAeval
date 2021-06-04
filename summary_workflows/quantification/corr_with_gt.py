import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bed', help='The BED file containing predictions merged with ground truth.')
args = parser.parse_args()

fname = args.bed
out = pd.read_csv(fname, delimiter='\t', header=None)

# initialize vectors for ground truth sites and prediction
vec_true = []
vec_pred = []

# multiple predicted sites for one ground truth?
multiple_predicted_sites = out.duplicated([6, 7, 8, 11], keep=False)

# iterate over matched sites
for (idx, row), is_mult in zip(out.iterrows(), multiple_predicted_sites):
    if is_mult:
        # sum up all prediction sites that were assigned to this ground truth site
        # needs to be implemented or can be skipped for now since these are usually only a few
        pass
    else:
        vec_true.append(row[10])
        # weighted expression in case there are multiple ground truth sites for one predicted site
        vec_pred.append(row[4]*row[12])

# correlation coefficient
r = pearsonr(vec_true, vec_pred)[0]

# export json
with open(fname[:-4] + '_corr.json', 'w') as f:
    json.dump({'correlation_coefficient': r}, f)

print('correlation coefficient: %.2f'%(r))
