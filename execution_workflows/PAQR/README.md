# PAQR
PAQR computes poly(A) site usage from RNA-seq data and a given reference poly(A) site set and is therefore only eligible for the **Quantification** benchmarking challenge.
## Rulegraph

![rulegraph](rulegraph.PAQR.png)

## Input & pre-processing
> **IMPORTANT:**   
PAQR requires samples, annotations and reference poly(A) sites files to all use the same chromosome naming scheme ("chr1" in gencode, "1" in ensembl) or else will fail **silently**!


1. Sample files  
`.bam` files of all samples have to be provided alongside their respective `.bai` index files (**in the same directory, with the same basename**)

2. Sample table  
A table specifying all samples to be analyzed has to be provided in `.tsv` format (NOT `.csv`) and contain the following columns:   

| ID | bam | condition |
| - | - | - |
| siRNA1 | path/to/bamfile1 | KO |
| siRNA2 | path/to/bamfile2 | KO |
| control1 | path/to/bamfile3 | CTRL |

>NOTE: If you already have a `.csv` file with the appropriate fields, you can convert it to `.tsv` BEFORE running the workflow with 
```
python ../../utils/csv2tsv/csv2tsv.py --csv config/samples.csv --tsv config/samples.tsv
```
 
See [here][sample-table] for an example input `samples.tsv` file. PAQR can be run on one or more samples, and one or more conditions. Conditions are only considered internally (for threshold calculations), the output is given for individual samples. PAQR does not support single-end data in reverse orientation.

This reference file has to be provided in `.bed` format with one poly(A) site per row. The site ID (column 4) has to be of the form `chr:site:strand` (e.g. "1:123456:+"), where "chr" is the chromosome, "site" is the representative site of the poly(A) site cluster, or the start coordinate in case of individual poly(A) sites, and "strand" is the strand on which the site is located. This format is based on [PolyASite][polyasite-web].    
Suitable reference poly(A) site files can be downloaded from [PolyASite][polyasite-web] for [human][human-pas], [mouse][mouse-pas] and [*C.elegans*][worm-pas]. The corresponding filename has to be specified in `config/config.PAQR.yaml` (see below).   
> **Note that PolyASite uses ensembl chromosome naming. For use with Gencode annotations, chromosome names have thus to be adjusted!**
```
# Terribly ugly but functional example for mouse

zcat atlas.clusters.2.0.GRCm38.96.bed.gz \
| awk 'BEGIN{OFS="\t"} {
    gsub(/^([0-9]+|[XY])/,"chr"$1,$1);
    gsub(/^MT/,"chrM",$1);
    split($4, id_fields, ":"); 
    gsub(/^([0-9]+|[XY])/,"chr"id_fields[1], id_fields[1]);
    gsub(/^MT/,"chrM", id_fields[1]);
    $4=id_fields[1]":"id_fields[2]":"id_fields[3];
    print}'\
| gzip -c > atlas.clusters.2.0.GRCm38.96.wchr.bed.gz

# Terribly ugly but functional example for human

zcat atlas.clusters.2.0.GRCh38.96.bed.gz \
| awk 'BEGIN{OFS="\t"} {
    gsub(/^([0-9]+|[XYM])/,"chr"$1,$1);
    split($4, id_fields, ":"); 
    gsub(/^([0-9]+|[XYM])/,"chr"id_fields[1], id_fields[1]);
    $4=id_fields[1]":"id_fields[2]":"id_fields[3];
    print}'\
| gzip -c > atlas.clusters.2.0.GRCh38.96.wchr.bed.gz
```

## Params

All parameters are specified (and explained) in `config/config.PAQR.yaml`. Most of the parameters don't have to be changed for a run with default behaviour, but do make sure the following ones are appropriate for your setup.

> NOTE: some parameters have to be specified more than once, with only slightly different names. Unfortunately, this cannot be avoided, as the workflow imports different individual modules, that all require to use the exact parameter names that are present in their respective published repositories. Re-wiring in order to decrease the overhead for the user has been performed as good as possible.
### Paths to input files
- `ref_PAS_file`
- `samples`
- `gtf`
- `PAQ_samples_table`


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
[human-pas]: <https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz>
[mouse-pas]: <(https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz)>
[worm-pas]: <https://polyasite.unibas.ch/download/atlas/2.0/WBcel235/atlas.clusters.2.0.WBcel235.bed.gz>
[sample-table]: config/samples.tsv
[ewf-readme-filenames]: ../README.md#output