## Matching PA sites to ground truth

The script `match_with_gt.py` uses bedtools window to assign ground truth PAS to the predictions. Each prediction is extended by `n` base pairs, depending on the parameter chosen for `window`.

Processing steps:
- run bedtools window for prediction file with itself
- if necessary, merge prediction sites that fall into the window, expression is summed
- run bedtools window for merged predictions with ground truth file
- find multiple ground truth sites overlapping one predicted site
	- weigths are added for the predicted site based on distance to the ground truth site
	- expression can be calculated by weight_i*expression_i for all prediction sites having the same ground truth site assigned
- find multiple predicted sites overlapping one ground truth site
	- weigths are added for the predicted sites based on distance to the ground truth site
	- expression can be calculated by sum(weight_i*expression_i) to get a single value matching the one ground truth site OR just summed without considering weight
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
python3 prediction.bed ground_truth.bed window

# example:
python3 match_with_gt.py adultCortex.PAPERCLIP.mm10.bed siControl_R1.MACEseq.mm10.bed 15
```

### Output:
- a copy of the prediction file `prediction_merged_<window>.bed` with rows merged, only in case there were overlapping sites given the window
- the output file containing ground truth matches `prediction_matched_<window>.bed`
