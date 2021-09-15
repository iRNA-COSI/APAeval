import subprocess
import pandas as pd
from shutil import which
from io import StringIO



def bedtools_window(f_GT, f_PD, window):
    """
    Python wrapper for bedtools window (by Leo S.)
    """
    
    # make sure bedtools can be called in current env
    assert which('bedtools') is not None, "bedtools not installed or not in PATH"

    # run bedtools window, capture output
    out = subprocess.run(['bedtools', 'window', '-sm',
                          '-w', str(window), 
                          '-a', f_GT, 
                          '-b', f_PD], 
                         capture_output=True, shell=False)
    assert out.returncode==0, "bedtools window run failed, check input files"

    # memory file-handle to pass output to pandas without writing to disk
    out_handle = StringIO(out.stdout.decode())
    
    # incase there were no sites returned (no overlap / all overlap in case of reverse=True)
    if not out.stdout.decode():
        out_df = pd.DataFrame()
    else:
        out_df = pd.read_csv(out_handle, delimiter='\t', header=None, dtype={0: str})
        
    out_df.rename({0: 'chrom_gt', 1: 'chromStart_gt', 2: 'chromEnd_gt', 3: 'name_gt', 4: 'score_gt', 5: 'strand_gt', 6: 'chrom_pd', 7: 'chromStart_pd', 8: 'chromEnd_pd', 9: 'name_pd', 10: 'score_pd', 11: 'strand_pd'}, axis=1, inplace=True)
   
    return out_df

def find_overlapping_sites(f_GT, f_PD, window):
    """
    Return sites in ground truth and test that are located withing the distance
    determined by window.
    """
    bt_window_out = bedtools_window(f_GT, f_PD, window)
    gt_true = list(set(bt_window_out["name_gt"]))
    pd_true = list(set(bt_window_out["name_pd"]))
    return gt_true, pd_true


def get_all_sites(bedfile):
    bed_df = pd.read_csv(bedfile, delimiter='\t', header=None, dtype={0: str})
    bed_df.rename({0: 'chrom', 1: 'chromStart', 2: 'chromEnd', 3: 'name', 4: 'score', 5: 'strand'}, axis=1, inplace=True)
    return list(bed_df["name"])

def count_sites(f_PD, f_GT, window):
    """
    Return TP (true positives), FP (false positives) and FN (false negatives).
    """
    
    gt_true, pd_true = find_overlapping_sites(f_GT, f_PD, window)
    
    gt_all_sites = get_all_sites(f_GT)
    pd_all_sites = get_all_sites(f_PD)
    
    gt_false = [x for x in gt_all_sites if x not in gt_true] # false negatives
    pd_false = [x for x in pd_all_sites if x not in pd_true] # false positives
    
    print(pd_true) # TODO: to be removed
    print(pd_false)
    print(gt_false)
    
    return len(pd_true), len(pd_false), len(gt_false)
 
def sensitivity(tp, fn):
    return tp / (tp + fn)

def fdr(tp, fp):
    return fp / (tp + fp)

def precision(tp, fp):
    return tp / (tp + fp)
    

   
def calculate_metrics(f_PD, f_GT, window):
    tp, fp, fn = count_sites(f_PD, f_GT, window)
    return sensitivity(tp, fn), fdr(tp, fp)
    