# PAQR
PAQR computes poly(A) site usage from RNA-seq data and a given reference poly(A) site set and is therefore only eligible for the **TPM & relative quantification** benchmarking challenges.
## Rulegraph

![rulegraph](rulegraph.PAQR.png)

## Input & pre-processing
> **IMPORTANT:**   
PAQR requires samples, annotations and reference poly(A) sites files to all use the same chromosome naming scheme ("chr1" in gencode, "1" in ensembl) or else will fail **silently**!


### 1. Sample files  
`.bam` files of all samples have to be provided alongside their respective `.bai` index files (**in the same directory, with the same basename**)

### 2. Sample table  
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

See [here][sample-table] for an example input `samples.tsv` file. PAQR can be run on one or more samples, and one or more conditions. Conditions are only considered internally (for threshold calculations), the output is given for individual samples. By default paired-end sequencing with read1 - reverse orientation, read2 - forward orientation is assumed. Single-stranded data with the reads in sense direction are processed properly too, but PAQR does not support single-end data in reverse orientation.

### 3. reference poly(A) site file   
This reference file has to be provided in `.bed` format with one poly(A) site per row. The site ID (column 4) has to be of the form `chr:site:strand` (e.g. "1:123456:+"), where "chr" is the chromosome, "site" is the representative site of the poly(A) site cluster, or the 1-based coordinate in case of individual poly(A) sites, and "strand" is the strand on which the site is located. This format is based on [PolyASite][polyasite-web].    
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

> NOTE: some parameters have to be specified more than once, with only slightly different names. Unfortunately, this cannot be avoided, as the workflow imports different individual modules, that all require to use the exact parameter names that are present in their respective published repositories. Re-wiring in order to decrease the overhead for the user has been performed as well as possible.
### Paths to input files
- `ref_PAS_file`: Path to reference poly(A) site file (format: `.bed` or `.bed.gz`)
- `tpas:gtf`: Path to genome annotation (format: `.gtf` or `.gtf.gz`)
- `paqr:PAQ_samples_table`: Path to samples table (format: `.tsv`)


### Data dependent
- `tpas:strandedness`: Specify whether the tandem PAS file should be created for use with stranded or unstranded data
- `paqr:PAQ_coverage_unstranded`: "no" for stranded data (has to match with `tpas:strandedness`)
- `paqr:PAQ_read_length`: avg read length of the samples

### Annotation dependent
- `biotype_key`: "transcript_biotype" if using ensembl annotations (also if "chr" has been prepended), "transcript_type" if using gencode


## Output & post-processing

APAeval relevant output: `tandem_pas_expression_normalized.tsv`, `tandem_pas_relative_usage.tsv` and `singular_pas_expression.tsv`. The latter is normalized in an additional rule outside the `PAQR module`, and subsequently both 'pas_expression' files are concatenated. The resulting `.tsv` file contains tpm for each sample from the samples table. This file is converted into the APAeval compatible 'format 02' bed file in a postprocessing rule. `tandem_pas_relative_usage.tsv` is directly converted into the APAeval compatible 'format 04' bed file in a post-processing rule. In both cases, single nucleotide PAS are reported for each PAS 'cluster' by parsing the representative position from the 'PolyASite ID'. See `method_workflow_file_specification.md` in the base `method_workflows` directory for a full specification.

The final output files are named `[SAMPLE]_[CHALLENGE_CODE]_[PARTICIPANT]_[OUTCODE].bed`, as specified in the ["method workflow README"][mwf-readme-filenames].

## Notes

The modules inside this workflow are loaded from the following repositories:

- [https://github.com/zavolanlab/tandem-pas](https://github.com/zavolanlab/tandem-pas)
- [https://github.com/zavolanlab/PAQR2](https://github.com/zavolanlab/PAQR2)

PAQR was developed to be compatible with specific downstream applications comparing distal and proximal PAS, and does thus natively only report expression of "tandem PAS" - PAS on terminal exons with at least two PAS. However, for APAeval also PAS from exons with only one PAS are relevant. Therefore the present workflow reports PAS expression from the PAQR "intermediate" files `singular_pas_expression.tsv` and `tandem_pas_expression_normalized.tsv`, instead of the original PAQR output `filtered_pas_expression.tsv`.


[polyasite-web]: <https://polyasite.unibas.ch/atlas>
[human-pas]: <https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz>
[mouse-pas]: <(https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz)>
[worm-pas]: <https://polyasite.unibas.ch/download/atlas/2.0/WBcel235/atlas.clusters.2.0.WBcel235.bed.gz>
[sample-table]: config/samples.tsv
[mwf-readme-filenames]: ../README.md#output
