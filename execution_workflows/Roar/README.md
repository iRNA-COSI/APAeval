# Roar
Roar detects the potential Alternative Polyadenylation (APA) events in the downstream 3'-UTR of a gene among two different biological conditions.

The Roar execution workflow is currently compatible with the following APAeval challenges:
- Differential challenge

(See [Notes](#Notes) section for clarification on exclusion from other challenges).


* Bioconductor: https://bioconductor.org/packages/release/bioc/html/roar.html
* Vignette: https://bioconductor.org/packages/release/bioc/vignettes/roar/inst/doc/roar.pdf
* Reference Manual: https://bioconductor.org/packages/release/bioc/manuals/roar/man/roar.pdf
* GitHub: https://github.com/vodkatad/roar


# Input files & parameters

The pipeline requires a 3 column, comma-separated sample table with the following headers:
- `sample_name` - unique name/identifier for the sample
- `bam` - path to BAM file of aligned reads for the given sample. It is expected (but not required) that an index file is present in the same location (suffixed with .bai)
- `condition` - string identifier for samples belonging to the same group

An example sample table can be found at `workflow/config/samples.csv` (works with APAeval test data).

The pipeline also requires additional reference files:
- **GTF file containing reference transcript models**
- **BED file containing reference poly(A) sites** (e.g. an atlas of polyA sites identified by 3'end-sequencing such as PolyASite, PolyADB or APASDB)
  - Poly(A) site coordinates should be defined as single nt length (see [Notes](#Notes) for explanation)

These reference files are used to generate Roar's custom APA site annotations. Since no public code is provided to generate these annotations, this is achieved using a custom script (`workflows/scripts/make_roar_annotations.py`) which follows instructions provided in the vignette. The script is currently only capable of generating '**[single APA type](https://github.com/vodkatad/roar/wiki/Identify-differential-APA-usage-from-RNA-seq-alignments)**' annotations.


The above files are provided via the pipeline's config file (e.g. `workflow/config/config.Roar.yaml`), which needs to be updated with run-specific information. The config files is currently set up to run directly with APAeval test data (see below for instructions).


# Running instructions

Activate the `apaeval_execution_workflows` conda environment ('environment' YAML file available at the base of the main repo at `apaeval_env.yaml`). If you haven't installed the environment, execute the following command (assuming you're in the same directory as this README):

```
conda env create -f ../../../apaeval_env.yaml`
```

Once installed, activate the environment with the command below:

```
conda activate apaeval_execution_workflows
```

Before running, you can perform a 'dry run' to check which steps will be run and where output files will be generated given the provided parameters and input files:

```
bash dryrun.sh
```

To run the workflow locally, you can use the provided wrapper script `run_local.sh`.

**Note: The run_local.sh script is currently set up to run with the APAeval test data**. If you have specified **absolute paths** in your sample sheet (e.g. `workflow/config/samples.csv`) or the config file (`workflow/config/config.Roar.yaml`), or have input data that is **not in the current directory**, you will need to modify Singularity bind arguments so the input files will be available to the container.

e.g. The path to the input GTF file is `/share/annotation/annotation.gtf`, and my current working directory is `/home/sam/DaPars2_snakemake/`. Modify the `--singularity-args` line in `run_local.sh` like below to ensure the file is available to the container:

```
--sigularity-args="--bind /share/" \
```

If you are satisfied with the bind arguments, you can run the workflow locally with the following command:

```
bash run_local.sh
```


# Output files

- Custom annotations generated for combination of provided GTF & BED file. In the main results directory this can be found at `roar_annotations.gtf`.
- TSV file storing results of Roar analysis (dataframe output by the *`totalResults()`* Roar function). In the main results directory this can be found at `roar_results.tsv`.
- TSV file compatible with the [differential usage challenge](https://github.com/iRNA-COSI/APAeval/blob/main/execution_workflows/execution_output_specification.md). In the main results directory this can be found named according to the `<CHALLENGECODE>_<PARTICIPANT>_<OUTCODE>.tsv` [specification](https://github.com/iRNA-COSI/APAeval/blob/main/execution_workflows/README.md#filenames) (each can be customised in the config file).


# Notes
> * Roar does not perform de-novo identification of polyA sites, instead relying on reference polyA sites. As such it is incapable of being compatible with the identification challenge
> * Roar does not output polyA site quantification as TPM values (only a relative change 'roar' (ratio of a ratio) metric and FPKM values), so is currently incompatible with the differential challenge
> * When defining the 3'end of the proximal polyA site bin (termed the 'PRE' region), the 'end' coordinate is lifted from the 'end' coordinate for the poly(A) site (*`_df_change_coordinate`* function in `workflow/scripts/make_roar_annotations.py`). If you use polyA site clusters (e.g. PolyASite), this will define the polyA site as the 3'end of the cluster, which may not be the most supported position by reads. The script will work if polyAsite coordinates are not single nucleotide, but note this may lead to small inaccuracies in the precise definition of the proximal polyA site.

# Citation

Grassi, E., Mariella, E., Lembo, A. et al. Roar: detecting alternative polyadenylation with standard mRNA sequencing libraries. BMC Bioinformatics 17, 423 (2016). https://doi.org/10.1186/s12859-016-1254-8

# Contact the Author
Elena Grassi - elena.grassi@unito.it, [GitHub](https://github.com/vodkatad)
