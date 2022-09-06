import subprocess
import os
from io import StringIO
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from shutil import which

#############################
# Genome/annotation related
#############################

def select_genome_file(file_name, genome_path):
    """Select the genome file according to the organism.

    Requires that the file_name contains an expression containing organism
    information, which will be matched against the genome_path directory.
    The format should be: name.mm10.ext or name.hg38extension.ext, with
    matching genome annotations: gencode.mm10.gtf and gencode.hg38extension.gtf.
    Note: no check for the extension (e.g. gtf) is done.

    Args:
        file_name (str): Name containing organism information. Supported: mm* and hg*.
        genome_path (str): directory containing genome annotations in gtf format.

    Returns:
        str with genome file path.
    """
    GENOME_STRINGS = ["mm", "hg"]
    SPLITSTRING = "."
    assert os.path.exists(genome_path), f"Genome annotation directory not found: {genome_path}"
    file_components =  file_name.split(SPLITSTRING)
    # search for genome
    for genome_string in GENOME_STRINGS:
        match = [comp for comp in file_components if genome_string in comp]
        if len(match) != 0:
            break
    if len(match) == 0:
        raise ValueError(f"No genome string: {GENOME_STRINGS} in file_name: {file_name} found.")
    # find all genome files in genome_path
    for f in os.listdir(genome_path):
        # find exact match in file
        genome_match = [f for comp in f.split(SPLITSTRING) if match[0] == comp]
        if len(genome_match) != 0:
            break
    if len(genome_match) == 0:
        raise ValueError(f"No genome string: {GENOME_STRINGS} in genome_path: {genome_path} found.")
    # return path to file
    return os.path.join(genome_path, genome_match[0])


def load_genome(genome_path):
    """Load genome annotation in gtf format.
    
    Requires feature 'gene', which is available in gencode.
    
    Args:
        genome_path: file path for gtf.
    
    Returns:
        pandas.df with genes.
    """
    assert os.path.exists(genome_path), f"Genome annotation not found: {genome_path}"
    # load gtf according to specifications (https://www.ensembl.org/info/website/upload/gff.html)
    genome = pd.read_csv(genome_path, sep="\t", header=None,
            comment="#",
            names=['seqname','source', 'feature','start','end','score','strand','frame','attribute'],
            usecols=['seqname','start','end','strand','feature'],
            dtype={'seqname': str, 'start': np.int64, 'end': np.int64, 'strand': str, 'feature': str})
    # Subset genome to only features gene
    subset = genome['feature'] == "gene"
    assert not subset.sum() == 0, f"No feature: 'gene' found in {genome_path}"
    genome.drop('feature', 1, inplace=True)
    return genome[subset]
    

def find_pas_genes(gene, df):
    """Find all PAS in given gene.

    Args:
        gene (dict): dict with keys 'seqname', 'start', 'end' and 'strand'.
        df (pandas.df): df with columns 'chrom_g' (str), 'chromStart_g' (int),
            'chromEnd_g' (int) and 'strand_g'.

    Returns:
        pd.Series: Array of indices indicating PAS for gene 'gene'.
    """
    subset = (df['chrom_g'] == gene['seqname']) & (df['strand_g'] == gene['strand'])
    subset = subset & (df['chromStart_g'] > gene['start']) & (df['chromEnd_g'] < gene['end'])
    return subset

#####################
# PAS matching
#####################


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
        out = pd.read_csv(out_handle, delimiter='\t', header=None, dtype={0: str,6: str})
        
    # label columns
    out.rename({0: 'chrom_p', 1: 'chromStart_p', 2: 'chromEnd_p', 3: 'name_p', 4: 'score_p', 5: 'strand_p', 6: 'chrom_g', 7: 'chromStart_g', 8: 'chromEnd_g', 9: 'name_g', 10: 'score_g', 11: 'strand_g'}, axis=1, inplace=True)
   
    return(out)

def find_weights(matched_sites_orig, pd_pos):
    """
    Calculate score weights for PD sites matched to multiple GT sites.
    If PD and GT sites are a perfect match, the weights are not split.
    If PD site between two GT sites, weights calculated based on distance.
    """
    matched_sites = matched_sites_orig.copy()
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

###############
# Metrics
###############


def precision(tp, fp):
    """
    Returns precision of prediction.
    """
    return tp / (tp + fp)

def sensitivity(tp, fn):
    """
    Returns sensitivity/TPR/recall of prediction.
    """
    return tp / (tp + fn)

def fdr(tp, fp):
    """
    Returns False Discovery Rate of prediction.
    """
    return fp / (tp + fp) 

def f1_score(precision, sensitivity):
    """
    Returns F1 score of prediction.
    """
    return  2 * (precision * sensitivity) / (precision + sensitivity)

def jaccard(tp, fp, fn):
    """
    Returns Jaccard index of prediction.
    """
    return  tp / (tp + fp + fn)

def corr_Pearson_with_gt(matched_sites):

    vec_true = matched_sites["score_g"]
    vec_pred = matched_sites["score_p"]

    # correlation coefficient
    r = pearsonr(vec_true, vec_pred)[0]

    return(r)

def corr_Spearman_with_gt(matched_sites):

    vec_true = matched_sites["score_g"]
    vec_pred = matched_sites["score_p"]

    # correlation coefficient
    r = spearmanr(vec_true, vec_pred)[0]

    return(r)

def relative_pas_usage(merged_bed_df, genome):
    """Compute correlation of relative PAS usage.

    1. find all PAS for a given gene, i.e. all ground truth PAS for a given gene or
        all predicted PAS that are matched to ground truth for a given gene.
    2. sum TPM values of all PAS
    3. calculate fraction for each PAS by dividing TPM_PAS by TPM_sum.

    Args:
        merged_bed_df (pandas.df): Table of matched prediction and ground truth PAS. From match_with_gt().
        genome (pandas.df): Table with gene positions.

    Returns:
        float: df of normalised and relative PAS values.
    """
    df = merged_bed_df.copy()
    # get list of PAS per gene
    pas_per_gene = [find_pas_genes(gene, df) for i, gene in genome.iterrows()]
    # remove PAS not covered by gene
    pas_per_gene = [pas for pas in pas_per_gene if len(np.where(pas)[0]) != 0]
    # compute fraction of PAS per gene
    normalised_dfs = [fraction_pas(gene_pas, df) for gene_pas in pas_per_gene]
    # concatenate list of pandas.df
    normalised_df = pd.concat(normalised_dfs, axis=0)
    return normalised_df

def fraction_pas(gene_pas, df):
    """Compute fraction for each PAS.
    Each PAS in the given gene is normalised by the sum of TPM values (column 'score'),
    separately for prediction and ground truth.

    Args:
        gene_pas (pandas.Series): bool vector indicating PAS for given gene.
        df (pandas.df): Table of matched prediction and ground truth PAS. From match_with_gt().
                        Columns 'score_p' and 'score_g' are used.
    
    Returns:
        pandas.df: 'df' for rows 'gene_pas' with 'score' columns as fractions.
    """
    pred_sum = df.loc[gene_pas, 'score_p'].sum()
    if pred_sum == 0:
        df.loc[gene_pas,'score_p'] = 0
    else:
        df.loc[gene_pas,'score_p'] = df.loc[gene_pas, 'score_p'] / pred_sum
    ground_sum = df.loc[gene_pas, 'score_g'].sum()
    if ground_sum == 0:
        df.loc[gene_pas,'score_g'] = 0
    else:
        df.loc[gene_pas,'score_g'] = df.loc[gene_pas,'score_g'] / ground_sum
    return df.loc[gene_pas,:]

