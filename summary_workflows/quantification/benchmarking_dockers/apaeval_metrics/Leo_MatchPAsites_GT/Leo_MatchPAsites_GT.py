import subprocess
from io import StringIO
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from shutil import which
import argparse

def bedtools_window(bed1, bed2, window, reverse=False):
    """
    Python wrapper for bedtools window.
    reverse: return only sites that have no match in the ground truth.
    """
    
    # make sure bedtools can be called in current env
    assert which('bedtools') is not None, "bedtools not installed or not in PATH"

    # run bedtools window, capture output
    if not reverse:
        out = subprocess.run(['bedtools', 'window', '-sm',
                              '-w', str(window), 
                              '-a', bed1, 
                              '-b', bed2], 
                             capture_output=True, shell=False)
    else:
        out = subprocess.run(['bedtools', 'window', '-sm', '-v',
                              '-w', str(window), 
                              '-a', bed1, 
                              '-b', bed2], 
                             capture_output=True, shell=False)

    assert out.returncode==0, "bedtools window run failed, check input files"

    # memory file-handle to pass output to pandas without writing to disk
    out_handle = StringIO(out.stdout.decode())
    
    # incase there were no sites returned (no overlap / all overlap in case of reverse=True)
    if not out.stdout.decode():
        out = pd.DataFrame()
    else:
        out = pd.read_csv(out_handle, delimiter='\t', header=None, dtype={0: str})
        
    # label columns
    out.rename({0: 'chrom_p', 1: 'chromStart_p', 2: 'chromEnd_p', 3: 'name_p', 4: 'score_p', 5: 'strand_p', 6: 'chrom_g', 7: 'chromStart_g', 8: 'chromEnd_g', 9: 'name_g', 10: 'score_g', 11: 'strand_g'}, axis=1, inplace=True)
   
    return(out)


def match_wtih_gt(f_PD, f_GT, window):

    # bedtools window with specified parameter
    out = bedtools_window(f_PD, f_GT, window)


    # initialize weight column with default 1 for a perfect match
    out['weight'] = [1]*len(out)

    # counter for logging
    weights_counter = 0

    # find rows in sample that are not unique 
    # -> extended predicted site overlaps multiple ground truth sites
    not_unq_mask = out.duplicated(['chrom_p', 'chromStart_p', 'chromEnd_p', 'strand_p'], keep=False)

    if np.sum(not_unq_mask) > 0: # otherwise there was no overlap
            
        not_unq = out[not_unq_mask]
        
        # remove non-unique sites
        out.drop(out.index[not_unq_mask], inplace=True, axis=0)
            
        # get the names of non-unique columns
        unq = np.unique([(i, j, k, l) for i, j, k, l in zip(not_unq['chrom_p'], not_unq['chromStart_p'], not_unq['chromEnd_p'], not_unq['strand_p'])], axis=0)

        for i in unq:
            ###### not implemented yet #######

            # split expression between ground truth sites, according to vicinity

            # re-name sites that have been split

            pass


    # find sites with no overlap given the window
    out_rev = bedtools_window(f_PD, f_GT, window, reverse=True)

    # for non-overlap sites, assign score 0 to ground truth and set other columns to values from prediction
    if not out_rev.empty:
        out_rev['chrom_g'], out_rev['chromStart_g'], out_rev['chromEnd_g'], out_rev['name_g'], out_rev['score_g'], out_rev['strand_g'], out_rev['weight'] = [out_rev['chrom_p'], out_rev['chromStart_p'], out_rev['chromEnd_p'], out_rev['name_p'], [0.0]*len(out_rev), out_rev['strand_p'], [1.0]*len(out_rev)]

    # add non-matched sites and matched sites together
    out = pd.concat([out, out_rev])

    # sort
    out.sort_values(by=['chrom_p', 'chromStart_p', 'chromEnd_p', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])

    # number of matched sites:
    n_matched_sites = len(out)-len(out_rev)
    # number of sites not matched:
    n_unmatched_sites = len(out_rev)
    
    return(out, n_matched_sites, n_unmatched_sites)

def corr_with_gt(out):

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

    return(r)


