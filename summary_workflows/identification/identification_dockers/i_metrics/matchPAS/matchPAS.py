import subprocess
import os
from io import StringIO
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from shutil import which

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

def find_weights(matched_sites, pd_pos):
    """
    Calculate score weights for PD sites matched to multiple GT sites.
    If PD and GT sites are a perfect match, the weights are not split.
    If PD site between two GT sites, weights calculated based on distance.
    """
    matched_sites.sort_values(by=["dist"], inplace=True)
    matched_sites = matched_sites.reset_index(drop=True)
    matched_sites["weight"] = 0.0

    # if the PD is not between GT sites, assume the closest GT is a perfect match
    if pd_pos < matched_sites["chromEnd_g"].min() or pd_pos > matched_sites["chromEnd_g"].max():
        matched_sites["dist"][matched_sites["dist"].idxmin()] = 0.0

    # calculate weight based on distance; a perfect match will have weight set to 1.0
    matched_sites["weight"][0:2] = 1 - (matched_sites["dist"][0:2] / np.sum(matched_sites["dist"][0:2])) # weight is equal to distance proportion
            
    return matched_sites

def split_pd_by_dist(matched_sites):
    """
    Identify PD sites matched to multiple GT sites and split the score between 2 sites based on distance
    """
    # find rows in sample that are not unique 
    # -> extended predicted site overlaps multiple ground truth sites
    
    unique_pd_cols = ['chrom_p', 'chromStart_p', 'chromEnd_p', 'strand_p']
    not_unq_mask = matched_sites.duplicated(unique_pd_cols, keep=False)
    

    if np.sum(not_unq_mask) > 0: # otherwise there was no overlap
            
        not_unq = matched_sites[not_unq_mask].copy()
       
        # remove non-unique sites
        matched_sites.drop(matched_sites.index[not_unq_mask], inplace=True, axis=0)
        
        # calculate distances between PD and GT sites
        not_unq["dist"] = abs(not_unq["chromEnd_p"] - not_unq["chromEnd_g"])

         # calculate weights based on the distance to two nearest GT sites
        not_unq = not_unq.groupby(unique_pd_cols).apply(lambda x: find_weights(x, x.name[2])).reset_index(drop=True)
        # apply weights to PD score
        not_unq["score_p"] = not_unq["score_p"] * not_unq["weight"]
        # rename sites that have been split
        not_unq["name_p"] = not_unq["name_p"] + "_" + not_unq["name_g"]

        # merge back to main df
        not_unq.drop(["dist", "weight"], axis=1, inplace=True)
        matched_sites = pd.concat([matched_sites, not_unq])

    return matched_sites

def merge_pd_by_gt(matched_sites):
    """
    Identify multiple PD sites matched to single GT site, merge and sum the scores.
    """
    
    # find GT rows in sample that are not unique - multiple PD matched to one GT
    # -> extended predicted site overlaps multiple ground truth sites
    unique_gt_cols = ['chrom_g', 'chromStart_g', 'chromEnd_g', 'strand_g']
    not_unq_GT_mask = matched_sites.duplicated(unique_gt_cols, keep=False)

    if np.sum(not_unq_GT_mask) > 0: # otherwise there was no overlap
        
        not_unq_GT = matched_sites[not_unq_GT_mask].copy()
        # remove non-unique sites
        matched_sites.drop(matched_sites.index[not_unq_GT_mask], inplace=True, axis=0)
        
        # sum scores for all PD sites matched to the same GT
        not_unq_GT["score_p"] = not_unq_GT.groupby(unique_gt_cols)["score_p"].transform("sum")
        # use only GT location for merged PD sites
        not_unq_GT[["chrom_p", "chromStart_p", "chromEnd_p"]] = not_unq_GT[["chrom_g", "chromStart_g", "chromEnd_g"]]
        not_unq_GT["name_p"] = not_unq_GT["name_g"] + "_merged"

        # drop duplicates and merge back to main df
        not_unq_GT.drop_duplicates(subset=["name_g"], keep="first", inplace=True)
        matched_sites = pd.concat([matched_sites, not_unq_GT])
    return matched_sites


def match_with_gt(f_PD, f_GT, window):

    # check for presence of participant input data
    assert os.path.exists(f_PD), "Participant file not found, please check input data."
    # check for presence of ground truth
    assert os.path.exists(f_GT), "Ground truth file not found, please check input data."
    
    # bedtools window with specified parameter
    out = bedtools_window(f_PD, f_GT, window)

    # splid PD sites that matched with multiple GT
    out = split_pd_by_dist(out)

    # find PD sites with no GT overlap given the window
    out_rev_PD = bedtools_window(f_PD, f_GT, window, reverse=True)

    # find GT sites with no PD overlap given the window, assign score 0 to ground truth and set other columns to values from prediction
    out_rev_GT = bedtools_window(f_GT, f_PD, window, reverse=True)
    if not out_rev_GT.empty:
        out_rev_GT[['chrom_g', 'chromStart_g', 'chromEnd_g', 'name_g', 'score_g', 'strand_g']] = out_rev_GT[['chrom_p', 'chromStart_p', 'chromEnd_p', 'name_p', 'score_p', 'strand_p']]
        out_rev_GT['score_p'] = [0.0]*len(out_rev_GT)
 
    # sort
    out.sort_values(by=['chrom_p', 'chromStart_p', 'chromEnd_p', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])

    # merge PD sites that matched with the same GT by summing their expression
    out = merge_pd_by_gt(out)

    # add GT sites with no PD match
    out_with_unmatched_GT = pd.concat([out, out_rev_GT])
    out_with_unmatched_GT.sort_values(by=['chrom_g', 'chromStart_g', 'chromEnd_g', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])

    # calculate total expression of non-matched PD sites
    if not out_rev_PD.empty:
        nonmatched_expression = (out_rev_PD["score_p"]).sum()
    else:
        nonmatched_expression = 0.0   
        
    return(out_with_unmatched_GT, nonmatched_expression)

def corr_with_gt(matched_sites):

    vec_true = matched_sites["score_g"]
    vec_pred = matched_sites["score_p"]

    # correlation coefficient
    r = pearsonr(vec_true, vec_pred)[0]

    return(r)
