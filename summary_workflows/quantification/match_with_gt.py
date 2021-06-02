import subprocess
from io import StringIO
import pandas as pd
import numpy as np
from shutil import which
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('a', help='The BED file containing predictions. MUST be BED6 format.')
parser.add_argument('b', help='The ground truth bed file. First 6 columns must be standard BED6, but can have additional columns appended.')
parser.add_argument('window', help='Number of bases to append to each side of the predicted site.', type=int)
#parser.add_argument('-o', help='output file directory') # not yet implemented!!
args = parser.parse_args()


f_GT = args.b
f_PD = args.a
window = args.window


def bedtools_window(bed1, bed2, window):
    """
    Python wrapper for bedtools window.
    """
    
    # make sure bedtools can be called in current env
    assert which('bedtools') is not None, "bedtools not installed or not in PATH"

    # run bedtools window, capture output
    out = subprocess.run(['bedtools', 'window', '-sm',
                          '-w', str(window), 
                          '-a', bed1, 
                          '-b', bed2], 
                         capture_output=True, shell=False)

    assert out.returncode==0, "bedtools window run failed, check input files"

    # memory file-handle to pass output to pandas without writing to disk
    out_handle = StringIO(out.stdout.decode())
    out = pd.read_csv(out_handle, delimiter='\t', header=None, dtype={0: str})
    
    return(out)

def overlap(coords1, coords2):
    """
    Calculates number of points overlapping when the input are two intervals as tuples.
    Input must be coords1=(int, int),  coords2=(int, int).
    """
    return(len(np.intersect1d(np.arange(coords1[0], coords1[1]), np.arange(coords2[0], coords2[1]))))


print('-------------------------------------------------')
## sum expression levels for predicted sites that overlap with each other given the window size
## assuming the prediction bed file has exactly 6 columns in the correct order, no header!

# find sites that overlap within the predict set itself when the window is extended
out = bedtools_window(f_PD, f_PD, window)

# for logging
len_dup = len(out)

not_unq_mask = out.duplicated([0, 1, 2, 5], keep=False)

MERGED=False
if np.sum(not_unq_mask) > 0: # otherwise there were no overlapping sites
    MERGED=True
    
    not_unq = out[not_unq_mask]

    # remove non-unique sites
    out.drop(out.index[not_unq_mask], inplace=True, axis=0)
    out.drop([6, 7, 8, 9, 10, 11], inplace=True, axis=1)

    # collapse non-unique sites, sum up expression, add back to bed file
    unq = np.unique([(i, j, k, l) for i, j, k, l in zip(not_unq[0], not_unq[1], not_unq[2], not_unq[5])], axis=0)

    for i in unq:
        non_unq_context = not_unq[(not_unq[0]==i[0]) & (not_unq[1]==int(i[1])) & (not_unq[2]==int(i[2])) & (not_unq[5]==i[3])]

        merged_start = min(list(non_unq_context[7]))
        merged_end = max(list(non_unq_context[8]))
        merged_exp = sum(non_unq_context[10])

        out = out.append({0: i[0], 1: merged_start, 2: merged_end, 3: list(non_unq_context[3])[0], 4: merged_exp, 5: i[3]}, ignore_index=True)

    # de-duplicate, sort
    out.drop_duplicates([0, 1, 2, 5], inplace=True)
    out.sort_values(by=[0, 1, 2], inplace=True, ascending=[True, True, True])
    
    f_PD_new = f_PD[:-4] + '_merged_' + str(window) + '.bed'
    out.to_csv(f_PD_new, sep='\t', header=None, index=False)
    
    print('Wrote intermediate file ' + f_PD_new + '.')
    print('Merged %i overlapping sites into %i sites.'%(len(unq), len_dup-len(unq)-len(out)))
    
else:
    print('No overlapping sites to merge within prediction file.')

print('-------------------------------------------------')


# matching ground truth with predicted (potentially merged) sites
## assuming the ground truth bed file has 6 columns, can have 1 additional column at the end with the gene name

if not MERGED:
    out = bedtools_window(f_PD, f_GT, window)
else:
    out = bedtools_window(f_PD_new, f_GT, window)

# initialize weight column with default 1 for a perfect match
out['weight'] = [1]*len(out)

# counter for logging
weights_counter = 0

# find rows in sample that are not unique 
# -> multiple ground truth sites overlapping one predicted site
not_unq_mask = out.duplicated([0, 1, 2, 5], keep=False)

if np.sum(not_unq_mask) > 0: # otherwise there was no overlap
        
    not_unq = out[not_unq_mask]
    
    # remove non-unique sites
    out.drop(out.index[not_unq_mask], inplace=True, axis=0)
        
    # collapse non-unique sites, sum up expression, add back to bed file
    unq = np.unique([(i, j, k, l) for i, j, k, l in zip(not_unq[0], not_unq[1], not_unq[2], not_unq[5])], axis=0)

    for i in unq:
        non_unq_context = not_unq[(not_unq[0]==i[0]) & (not_unq[1]==int(i[1])) & (not_unq[2]==int(i[2])) & (not_unq[5]==i[3])].copy()
        
        olps = []
        for idx, gt_site in non_unq_context.iterrows():
            frac_overlap = overlap((gt_site[1]-window, gt_site[2]+window), (gt_site[7]-window, gt_site[8]+window))/(gt_site[2]-gt_site[1]+2*window)
            olps.append(frac_overlap)
            weights_counter += 1
        
        weights = np.array(olps)/np.sum(olps)
        non_unq_context['weight'] = weights
        
        # add back the removed sites including weights
        out = out.append(non_unq_context)

    # sort
    out.sort_values(by=[0, 1, 2, 7], inplace=True, ascending=[True, True, True, True])
    
# find rows in ground truth that are not unique 
# -> multiple predicted sites overlapping one ground truth site
not_unq_mask = out.duplicated([6, 7, 8, 11, 'weight'], keep=False)

if np.sum(not_unq_mask) > 0: # otherwise there was no overlap
        
    not_unq = out[not_unq_mask]
    
    # remove non-unique sites
    out.drop(out.index[not_unq_mask], inplace=True, axis=0)

    # collapse non-unique sites, sum up expression, add back to bed file
    unq = np.unique([(i, j, k, l) for i, j, k, l in zip(not_unq[6], not_unq[7], not_unq[8], not_unq[11])], axis=0)
    
    for i in unq:
        non_unq_context = not_unq[(not_unq[6]==i[0]) & (not_unq[7]==int(i[1])) & (not_unq[8]==int(i[2])) & (not_unq[11]==i[3])].copy()
        
        olps = []
        for idx, gt_site in non_unq_context.iterrows():
            frac_overlap = overlap((gt_site[1]-window, gt_site[2]+window), (gt_site[7]-window, gt_site[8]+window))/(gt_site[2]-gt_site[1]+2*window)
            olps.append(frac_overlap)
            weights_counter += 1
        
        weights = np.array(olps)/np.sum(olps)
        non_unq_context['weight'] = weights
        
        # add back the removed sites including weights
        out = out.append(non_unq_context)

    # sort
    out.sort_values(by=[0, 1, 2, 7], inplace=True, ascending=[True, True, True, True])

f_PD_matched = f_PD[:-4] + '_matched_' + str(window) + '.bed'
out.to_csv(f_PD_matched, sep='\t', header=None, index=False)
print('Wrote output to file ' + f_PD_matched)
print('Matched %i sites in total.'%(len(out)))
print('Added weights to %i overlapping sites.'%(weights_counter))

print('-------------------------------------------------')
