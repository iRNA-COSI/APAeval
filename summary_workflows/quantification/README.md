## Matching PA sites to ground truth

The script `match_with_gt.py` uses bedtools window to assign ground truth PAS to the predictions. Each prediction is extended by `n` base pairs in both directions, depending on the parameter chosen for `window_size`.

The script `corr_with_gt.py` calculates the correlation coefficient between prediction and matched ground truth quantification values. The input must be a BED file containing the columns below. Usage: `python3 corr_with_gt.py -o <output.json> <prediction_matched.bed>`

Processing steps:
- run bedtools window for predictions and ground truth file
- run bedtools window for predictions and ground truth file with -v option to retrieve PAS without a match
- find multiple ground truth sites overlapping one predicted site (currently these cases are excluded)
        - predicted sites are split and the expression weighted based on the overlap with the ground truth sites
- write output file with the following columns:
        1. chrom prediction
        2. chromStart prediction
        3. chromEnd prediction
        4. name prediction
        5. score prediction (expression value)
        6. strand prediction
        7. chrom ground truth
        8. chromStart ground truth
        9. chromEnd ground truth
        10. name ground truth
        11. score ground truth (additional columns go after this one, such as ground truth gene_ID)
        12. strand ground truth
        13. weight of prediction expression

### Usage:

```bash
python3 match_with_gt.py -o output.bed prediction.bed ground_truth.bed window_size

# example:
python3 match_with_gt.py adultCortex.PAPERCLIP.mm10.bed siControl_R1.MACEseq.mm10.bed 15
```

### Output:
- contains all rows of the prediction file, with cases where a predicted site matches multiple ground truth sites removed
- when no overlapping ground truth is found, the ground truth columns are transferred from the prediction and the expression is set to 0

### To do:
- handle cases where multiple ground truth sites are found overlapping one predicted site
