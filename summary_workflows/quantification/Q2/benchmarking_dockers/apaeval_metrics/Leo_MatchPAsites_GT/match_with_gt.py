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


def bedtools_window(bed1, bed2, window, reverse=False):
    """
    Python wrapper for bedtools window.
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
    out = pd.read_csv(out_handle, delimiter='\t', header=None, dtype={0: str})
    
    return(out)

def overlap(coords1, coords2):
    """
    Calculates number of points overlapping when the input are two intervals as tuples.
    Input must be coords1=(int, int),  coords2=(int, int).
    """
    return(len(np.intersect1d(np.arange(coords1[0], coords1[1]), np.arange(coords2[0], coords2[1]))))


print('-------------------------------------------------')


# bedtools window with specified parameter
out = bedtools_window(f_PD, f_GT, window)
# label columns
out.rename({0: 'chrom_p', 1: 'chromStart_p', 2: 'chromEnd_p', 3: 'name_p', 4: 'score_p', 5: 'strand_p', 6: 'chrom_g', 7: 'chromStart_g', 8: 'chromEnd_g', 9: 'name_g', 10: 'score_g', 11: 'strand_g'}, axis=1, inplace=True)


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
        # split expression between ground truth sites, according to vicinity

        # re-name sites that have been split


    
        pass



        # non_unq_context = not_unq[(not_unq['chrom_p']==i[0]) & (not_unq['chromStart_p']==int(i[1])) & (not_unq['chromEnd_p']==int(i[2])) & (not_unq['strand_p']==i[3])].copy()
        
        # olps = []
        # for idx, gt_site in non_unq_context.iterrows():
        #     frac_overlap = overlap((gt_site['chromStart_p']-window, gt_site['chromEnd_p']+window), (gt_site['chromStart_g']-window, gt_site['chromEnd_g']+window))/(gt_site['chromEnd_p']-gt_site['chromStart_p']+2*window)
        #     olps.append(frac_overlap)
        #     weights_counter += 1
        
        # weights = np.array(olps)/np.sum(olps)
        # non_unq_context['weight'] = weights
        
        # # add back the removed sites including weights
        # out = out.append(non_unq_context)

# find sites with no overlap given the window
out_rev = bedtools_window(f_PD, f_GT, window, reverse=True)
out_rev.rename({0: 'chrom_p', 1: 'chromStart_p', 2: 'chromEnd_p', 3: 'name_p', 4: 'score_p', 5: 'strand_p', 6: 'chrom_g', 7: 'chromStart_g', 8: 'chromEnd_g', 9: 'name_g', 10: 'score_g', 11: 'strand_g'}, axis=1, inplace=True)

# for non-overlap sites, assign score 0 to ground truth and set other columns to values from prediction
out_rev['chrom_g'], out_rev['chromStart_g'], out_rev['chromEnd_g'], out_rev['name_g'], out_rev['score_g'], out_rev['strand_g'], out_rev['weight'] = [out_rev['chrom_p'], out_rev['chromStart_p'], out_rev['chromEnd_p'], out_rev['name_p'], [0.0]*len(out_rev), out_rev['strand_p'], [1.0]*len(out_rev)]

# add non-matched sites and matched sites together
out = pd.concat([out, out_rev])

# sort
out.sort_values(by=['chrom_p', 'chromStart_p', 'chromEnd_p', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])
    
f_PD_matched = f_PD[:-4] + '_matched_' + str(window) + '.bed'
out.to_csv(f_PD_matched, sep='\t', header=None, index=False)
print('Wrote output to file ' + f_PD_matched)
print('Matched %i sites in total.'%(len(out)-len(out_rev)))
print('Added weights to %i overlapping sites.'%(weights_counter))

print('-------------------------------------------------')
