# PAQR
## Rulegraph

![rulegraph](rulegraph.PAQR.png)

## Input & pre-processing
> IMPORTANT:   
PAQR requires samples, annotations and reference poly(A) sites files to all use the same chromosome naming scheme ("chr1" in gencode, "1" in ensembl) or else will fail.


1. Sample files  
`.bam` files of all samples have to be provided alongside their respective `.bai` index files (**in the same directory, with the same basename**)

2. Sample table  
A table specifying all samples to be analyzed has to be provided in `.tsv` format (not `.csv`) and contain the following columns:   

| ID | bam | condition |
| - | - | - |
| siRNA1 | path/to/bamfile1 | KO |
| siRNA2 | path/to/bamfile2 | KO |
| control1 | path/to/bamfile3 | CTRL |
 
See [here][sample-table] for an example. PAQR can be run on one or more samples, and one or more conditions. Conditions are only considered internally (for threshold calculations), the output is given for individual samples.

3. reference poly(A) site file   
This file has to be provided in `.bed` format with one poly(A) site per row. The site ID (column 4) has to be of the form `chr:site:strand` (e.g. "1:123456:+"), where "chr" is the chromosome, "site" is the representative site of the poly(A) site cluster, or the start coordinate in case of individual poly(A) sites, and "strand" is the strand on which the site is located. This format is based on [PolyASite][polyasite-web].
## Params

All parameters are specified (and explained) in `config/config.yaml`. Most of the parameters don't have to be changed for a run with default behaviour, but do make sure the following ones are appropriate for your setup.

> NOTE: some parameters have to be specified more than once, with only slightly different names. Unfortunately, this cannot be avoided, as the workflow imports different individual modules, that all require to use the exact parameter names that are present in their respective published repositories.
### Paths to input files
- `samples`
- `polyasite`
- `gtf`
- `PAQ_samples_table`
- `PAQ_tandem_pas`: This file is created as output of the tandem PAS module and its location and name follows the format `[outdir]/[atlas_version]_tandem_pas.terminal_exons.[strandedness].bed`, where expressions in brackets correspond to the respective values in `config/config.yaml`.

### Data dependent
- `strandedness`: Specify whether the tandem PAS file should be created for use with stranded or unstranded data
- `PAQ_coverage_unstranded`: "no" for stranded data
- `PAQ_read_length`: avg read length of the samples

### Annotation dependent
- `biotype_key`: "transcript_biotype" if using ensembl annotations, "transcript_type" if using gencode


## Output & post-processing

APAeval relevant output: `filtered_pas_expression.tsv`, which contains tpm for each sample from the samples table. This file is converted into the APAeval compatible bed format in a postprocessing rule.   

The final output files are named `[SAMPLE]_[CHALLENGE_CODE]_[PARTICIPANT]_[OUTCODE].bed`, as specified in the ["execution workflow README"][ewf-readme-filenames].

## Notes

The modules inside this workflow are loaded from the following repositories:

- [https://github.com/zavolanlab/tandem-pas](https://github.com/zavolanlab/tandem-pas)
- [https://github.com/zavolanlab/PAQR2](https://github.com/zavolanlab/PAQR2)



[polyasite-web]: <https://polyasite.unibas.ch/atlas>
[sample-table]: config/samples.tsv
[ewf-readme-filenames]: ../README.md#output