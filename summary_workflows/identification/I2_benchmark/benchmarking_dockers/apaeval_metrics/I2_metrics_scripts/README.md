## Functions for calculating sensitivity and FDR for I2 benchmark

The functions in `i2_metrics.py` use `bedtools window` for matching PAS from predictions to ground truth datasets for a given value of `window` parameter.

Any predicted PAS that fall within `window` from ground-truth PAS are assumed to be **true positives** even if:

- there are multiple predicted PAS for one ground truth PAS
- predicted PAS is matched to more than one ground truth PAS (possible for large window size)

**False positives** are all predicted PAS that don't match with ground-truth PAS for given window size.

**False negatives** are all ground truth PAS that don't have matching predicted PAS for a given window size.

Based on these values, the following metrics can be calculated:

1. Sensitivity = (TP/(TP+FN))
2. FDR = (FP/(TP+FP))

The scripts assume the following BED file fields for both predicted and ground-truth PAS:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is the same as starting position
- **name** - defines the **unique** name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-"

**Note**: for now the functions use the **name** field, and the assumption is that there is a unique name for each site in both datasets.




