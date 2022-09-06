# APAeval comparison of PAS between predictions and ground truths

This module contains functions for
- matching of poly(A) sites between prediction (PD) and ground truth (GT)
- calculation of performance metrics for benchmarking the prediction tools

Below is an overview of the contents of the module. If available, the documentation for individual functions can be found [here][module-md].



## PAS matching
- [Bedtools window][bedtools-window] is used to match poly(A) sites between GT and PD, allowing an imprecision of the size of the specified window upstream and downstream of the GT site.
- If a PD site matches multiple GT sites its score is split between the GT sites according to distance:
  - weights are added for the PD site based on distance to the GT site if the predicted site is between two GT sites; if it is outside, it is assumed that the closest GT is a perfect match
  - for perfect matches the weight is set to 1, which means that the other overlapping sites get weight 0
  - expression can be calculated by weight_i*expression_i for all PD sites having the same GT site assigned
- If multiple PD sites match one GT they are merged and their score is summed without applying distance based weights
- False positives are PD sites without a GT match determined by [`bedtools_window -v`][bedtools-window]; expression of GT is set to 0 for those sites
- False negatives are GT sites without a PD match determined by [`bedtools_window -v`][bedtools-window]; expression of PD is set to 0 for those sites
- matched sites, FP and FN are returned as three separate pandas dataframes with the following columns each:
	1. chromosome PD
	2. start PD
	3. end PD
	4. name PD
	5. expression PD
	6. strand PD
	7. chromosome GT
	8. start GT
	9. end GT
	10. name GT
	11. expression GT 
	12. strand GT

> ATTENTION: If a tool reports all PAS present in its reference DB regardless of expression, the *number* of FP for that tool will be overrated (whereas the *expression* of FP remains correct). For identification evaluations filtering for sites with non-zero expression thus has to be performed.

## Metrics calculations
- Pearson correlation coefficient
- Spearman correlation coefficient
- Correlation of relative PAS usage per gene

## Other functions
- load annotation of genes
- find PAS that belong to given gene
- normalized expression of PAS in gene


## Developer notes
- Locally build the module for debugging: `pip install -e .` into the active APAeval conda environment.
- For creating the markdown documentation of the module: 
	```bash
	# inside activated APAeval conda env
	pip install pdoc3

	# cd into utils/apaeval/
	pdoc -o . src/apaeval/main.py
	```



## Origin
This module is based on original APAeval2021 hackathon work done by [Leo Schärfen](https://github.com/lschaerfen):

- [`match_with_gt.py`](https://github.com/iRNA-COSI/APAeval/blob/9a17c11dd6239969feb092d687ac7e206043c8d6/summary_workflows/quantification/match_with_gt.py)
- [`corr_with_gt.py`](https://github.com/iRNA-COSI/APAeval/blob/9a17c11dd6239969feb092d687ac7e206043c8d6/summary_workflows/quantification/corr_with_gt.py)

[//]: # (References)
[module-md]: ./main.md
[bedtools-window]: https://bedtools.readthedocs.io/en/latest/content/tools/window.html