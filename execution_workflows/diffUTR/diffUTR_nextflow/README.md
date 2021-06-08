# diffUTR
The diffUTR R package streamlines differential exon usage (DEU) analyses, and leverages existing DEU tools and alternative poly-adenylation site databases to enable differential 3' UTR usage analysis.


#### Links:
- Paper: [Streamlining differential exon and 3′ UTR usage with diffUTR](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04114-7)
- [Vignette](https://bioconductor.org/packages/release/bioc/vignettes/diffUTR/inst/doc/diffUTR.html)
- GitHub [repo](https://github.com/ETHZ-INS/diffUTR)

## Input & pre-processing

Sample table: view [here](input.csv)

Specified in the `input.csv`:

- Standard GTF (e.g. GENCODE, no __"chr"__ in the chrom column)
- PolyA site annotation (e.g. a bed file from PolyASite, PolyA_DB, can be downloaded inside the package) [Not in the `test_data/` dir yet]
- BAM files
- Method specification (specify `DEXSeq`, `limma`, or `edgeR`):
  - DEXSeq
  - limma's `diffSplice` (by default, as it is the one used in their vignette)
  - edgeR's `diffSpliceDGE`
- Sequencing information: 
  - strandSpecific: specify `unstranded`, `stranded`, or `reversely stranded`
  - single or paired end: specify `SE` or `PE`
- Analysis type:
  - Differential exon
  - Differential UTR


## Params

Descriptions see above. Example see the [`input.csv`](input.csv) file

```
--input input.csv
```


## Output & post-processing

diffUTR output

## Notes

- The samples need to have replicates (cannot be 1v1)
- The "chr" in the chromosome column in GTF has to be removed, if not there will be an error.
- There is a bug, technical issue, in this package, I opened a [PR](https://github.com/ETHZ-INS/diffUTR/pull/4) in the repo
- The covariate used in the differential analysis is fixed here, i.e. `condition` (e.g. limma's linear model: `design ~ condition`)
