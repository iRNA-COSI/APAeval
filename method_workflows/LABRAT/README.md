# Method workflow for LABRAT

`LABRAT` is for quantifying the alternative usage of polyadenylation and cleavage sites in RNAseq data. You can find the [publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07781-1).
source code: https://github.com/TaliaferroLab/LABRAT

The quantification step of LABRAT does not qualify for the benchmark quantification, which requires TPM values for each polyA site. LABRAT only provides TPM per transcript. 
However, the output could be generated in the future by collecting the polyA site coordinates. 

## Input & pre-processing
- samplesheet.csv with the following columns (more_samplesheets/samplesheet_paired.csv was used for testing locally)
  - sample -- unique names of the samples
  - fastq1 -- single-end fastq or paired-end fastq1
  - fastq2 -- paired-end fastq2 (if single-end, please leave it empty)
  - gff -- gencode gff file (tested Ensembl does not work due to the lack of the gene_type field)
  - fasta -- genome fasta file 
  - condition -- condition of the samples (requires at least two conditions)

## Params

- `--input` samplesheet.csv
- `--run_quantification` for running quantification
- `--run_differential` for running differential analysis
- `--labrat_quant_dir` location of the quantification results directory. Please provide this if `--run_differential` but not `--run_quantification`
- `--conditionA` the first condition in the condition column of the samplesheet.csv. Please provide this if `--run_differential`
- `--conditionB` the second condition in the condition column of the samplesheet.csv. Please provide this if `--run_differential`

The Nextflow workflow contains two profiles: *docker* and *singularity*. Either can be used for exeuction and they are defined in `nextflow.config`.

## example commands
#### running quantification only
```
nextflow main.nf -profile docker --input <samplesheet.csv> --run_quantification
```
#### running differential analysis only
```
nextflow main.nf -profile docker --input <samplesheet.csv> --run_differential --conditionA <first condition> --conditionB <second condition> --labrat_quant_dir <full_path to labrat_quantification_results_dir>
```
#### running both quantification and differential analysis
```
nextflow main.nf -profile docker --input <samplesheet.csv> --run_quantification --run_differential --conditionA <first_condition> --conditionB <second_condition>
```

#### Run test data

The workflow does not contain the setup to run the test data in `tests/test_data` as the genome sequence is necessary. 
The following steps are meant to setup the test data in order to be comparable with other method workflows.

   * Create `siControl_R1_2genes.fastq.gz` and `siSrsf3_R1_2genes.fastq.gz` from corresponding BAM files with `samtools bam2fq` as described here: `tests/test_data/README.md`.
   * Use `gencode_2genes_Chr_prefix.vM18.annotation.gff3` as genome annotation (gff).
   * Download genome sequence of *mm10* from Gencode, preferentially only *chr6* and *chr16* as only those are used in the test.
   * Provide the files with absolute paths in a sample sheet, by e.g. adjusting `more_samplesheets/samplesheet_single.csv`.
      * Leave column *fastq2* empty
      * sample: *siControl_R1* and *siSrsf3_R1*
      * condition: *siControl_* and *siSrsf3*
   * Run the quantification and differential analysis with (use actual sample sheet):
   ```
   nextflow main.nf --profile docker --input more_samplesheets/samplesheet_single.csv \
      --run_quantification --run_differential \
      --conditionA siControl --conditionB siSrsf3
   ```

## Output & post-processing

`LABRAT` consists of the following three steps:
- `makeTFfasta` -- makes a fasta file of transcripts that will later be quantified by salmon
- `runSalmon` -- quantifies transcript abudance with salmon (please see Notes below)
   - each sample will have its own output directory
   - `quant.sf` is used for generating the quantification bed file but without `chromStart` and `chromEnd` inforamtion (the starting and ending position of the feature in the chromosome.)
- `calculatepsi` -- calculates the 𝜓 values for each gene in each sample and the identify genes that show significantly different 𝜓 values across conditions.
   - The output `LABRAT.psis.pval` is used for generation the differential usage tsv file

## Notes

The `LABRAT.py --mode runSalmon` step is replaced by running `salmon` directly with the [`salmon index`](https://github.com/TaliaferroLab/LABRAT/blob/73b9ddd6b7922a349419b49b2ed25e93e9f261f0/LABRAT.py#L877) and [`salmon quant` commands indicated in `LABRAT.py`](https://github.com/TaliaferroLab/LABRAT/blob/73b9ddd6b7922a349419b49b2ed25e93e9f261f0/LABRAT.py#L354) due to c++ compling problems of salmon in LABRAT's BioContainer.
