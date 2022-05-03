## Matching PA sites to ground truth

The script `matchPAS.py` contains the metric calculating functions that were in the following scripts [Leo Schärfen](https://github.com/lschaerfen) wrote: 

- [`match_with_gt.py`](https://github.com/iRNA-COSI/APAeval/blob/9a17c11dd6239969feb092d687ac7e206043c8d6/summary_workflows/quantification/match_with_gt.py) uses bedtools window to assign ground truth PAS to the predictions. Each prediction is extended by `n` base pairs, depending on the parameter chosen for `window`.
- [`corr_with_gt.py`](https://github.com/iRNA-COSI/APAeval/blob/9a17c11dd6239969feb092d687ac7e206043c8d6/summary_workflows/quantification/corr_with_gt.py) calculates the correlation coefficient between prediction and matched ground truth quantification values. The input must be a BED file containing the columns below. Usage: `python3 corr_with_gt.py <prediction_merged.bed>`

Processing steps:
- run bedtools window for prediction file with itself
- if necessary, merge prediction sites that fall into the window, expression is summed
- run bedtools window for merged predictions with ground truth file
- find multiple ground truth sites overlapping one predicted site
	- weights are added for the predicted site based on distance to the ground truth site if the predicted site is between two ground truth sites; if the predicted site matches one of the ground truth sites perfectly, the weight for this site is equal to 1 and for other sites it is equal to 0
	- expression can be calculated by weight_i*expression_i for all prediction sites having the same ground truth site assigned
- find multiple predicted sites overlapping one ground truth site
	- weights are added for the predicted sites based on distance to the ground truth site
	- expression can be calculated by sum(weight_i*expression_i) to get a single value matching the one ground truth site OR just summed without considering weight
- ground truth sites that do not have any matching predicted sites have expression for the predicted site set to 0 and are used for correlation calculations
- predicted sites that do not have any matching ground truth sites are handled separately
- write output file with the following columns:
	1. chromosome prediction
	2. start prediction
	3. end prediction
	4. name prediction
	5. expression prediction
	6. strand prediction
	7. chromosome ground truth
	8. start ground truth
	9. end ground truth
	10. name ground truth
	11. expression ground truth (additional columns go after this one, such as ground truth gene_ID)
	12. weight of prediction expression

### Usage:

```bash
python3 matchPAS.py prediction.bed ground_truth.bed window

# example:
python3 matchPAS.py adultCortex.PAPERCLIP.mm10.bed siControl_R1.MACEseq.mm10.bed 15
```

### Output:
- a copy of the prediction file `prediction_merged_<window>.bed` with rows merged, only in case there were overlapping sites given the window
- the output file containing ground truth matches `prediction_matched_<window>.bed`


### Issues:
- currently, predicted sites that have no overlapping ground truth site within the window are discarded
- there is no option in `bedtools window` to keep those, except for `-v` which requires another  `bedtools window` run
- this needs to be addressed because it would favor tools that find less sites but with higher accuracy
