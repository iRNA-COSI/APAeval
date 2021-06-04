# {Snakemake pipeline template}

## How to create your snakemake pipeline for APAeval
Read this section, but do NOT include it your final README.
> * Copy the whole directory `snakemake` with its' subfolders to suitable location.
> * Rename `snakemake` to `[METHOD]_snakemake`, [METHOD] being the name of the tool you're building a workflow for.
> * Adjust names and contents of files that appear in this [NOTATION].
> * Place any scripts or subworkflows you're going to use into the directories `workflow/scripts` and `workflow/rules`, respectively.
> * Write your workflow in `Snakefile`.
> * Here are the [snakemake docs](https://snakemake.readthedocs.io/en/stable/index.html)
> * Create the [conda](https://docs.conda.io/en/latest/) environment named `snakemake` with `conda env create -f snakemake.yaml`. Alternatively, create an own virtual environment, but ensure to use same versions as in `snakemake.yaml`.
> * General good practices:
>     * Test your code with [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#best-practices).
>     * One (shell) command per rule. 
>     * If samples differ in a meaningful way (e.g. single end and paired end samples), it might be better to write subworkflows within `workflow/rules`. 
>     * If applicable, each rule has it's own [Docker container](https://www.docker.com/resources/what-container) or conda environment.
> * There are some shell scripts which can be used to start a snakemake run. Adjust the name of the `config file` in these scripts:
>     * `dryrun.sh`
>     * `rulegraph.sh` (Also adjust name of output `.png`)
>     * `run_local.sh`
> * Check out the pilot benchmark at `tests/pilot_benchmark/snakemake` for a running example. It illustrates how the execution workflow can look in practice. 
> * Adjust this `README.md`: Delete this description section and populate the sections below with your awesome experience ;)

# DaPars2 {Method name as specified in Algorithm table.}

## Rulegraph

{The rulegraph gives an overview of the steps of the workflow. To obtain it, adapt and run the `rulegraph.sh` script!}

![rulegraph](rulegraph.[METHOD].png)

## Input & pre-processing
> * Note: the chromosome name conatin no 'chr' prefix in the files below, this workflow will add chr to each file.
>   * Bam files
>   * Bed files (sample)
>   * Gtf annotation
>   * Gff3 annotation
>   * Bed file of the genome

{Describe input files and how they will be processed in order for the method to work. Describe how sample tables have to look like, and any other input that is needed (e.g. genome).}

## Params
>   * Coverage Threshold
{Describe parameters needed by your METHOD.}

## Output & post-processing
> * Output directory chr1-22, chrX, chrY for each sample. Each directory contains DaPars2 results.
> * 01_DaPars2.bed

{Describe output files and postprocessing steps if necessary.}

## Notes
> * In `DaPars2_Multi_Sample_Multi_Chr.py` mode, DaPars2 only processes samples with chromosome name with 'chr' prefix. However, if running with chromosome specific analysis, `Dapars2_Multi_Sample.py` mode, one can assign the name in the input samples. This workflow executes `DaPars2_Multi_Sample_Multi_Chr.py` only.
> * In [DaPars2 documentation](http://bioinfo.szbl.ac.cn/DaPars2/DaPars2.html), the input files are in wiggle format; however, the testing data that it provides is in bedgraph formant. 
{Notes about the METHOD. 
e.g. Did you have to adjust the method's soure code?
}
