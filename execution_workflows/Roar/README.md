# Roar
Roar detects the potential Alternative Polyadenylation (APA) events in the downstream 3'-UTR of a gene among two different biological conditions.

# Installation
Roar can be installed from GitHub & BiocManager

* GitHub: https://github.com/vodkatad/roar
* BiocBioconductor: https://bioconductor.org/packages/release/bioc/html/roar.html

# User manual
Documentation for Roar is available at:

* Tutorial: https://bioconductor.org/packages/release/bioc/vignettes/roar/inst/doc/roar.pdf
* User Manual: https://bioconductor.org/packages/release/bioc/manuals/roar/man/roar.pdf

# Environments

**A Docker image still needs to be made for this**. All the needed packages and tools needed to run this are present in the conda `.yaml`, except for Roar itself which can be installed from GitHub/BioConductor (above).

# Workflow (Snakemake)

Workflow is written in **Snakemake**, saved as a `workflow/Snakefile`. There is just one rule **run_roar_rscript**, which runs the R script `workflows/scripts/Roar.R` which does all the pre-processing, runs roar (multiple APA analysis, see the roar tutorial section 3.6 for more detauls) and preprocesses everything. It takes as input the R script file, a `.gtf` annotation file of 3' UTR APA sites for your organism of choice. The `.gtf` files are found in `gtfs` and were downloaded from https://github.com/vodkatad/roar/wiki/Identify-differential-APA-usage-from-RNA-seq-alignments and are the ones for the multiple APA analysis, generated from either PolyADB or APAsdb. It also takes at 2 pairs of `.bam` files of 2 different levels, located in folders that are named after levels, the names of which are present in the "sample_type" column of the `config/samples.csv` file. For example with the test data, they'd be present in `srsf3` and `control`.


# Output files

The workflow outputs 2 formats (https://github.com/iRNA-COSI/APAeval/blob/main/execution_workflows/execution_output_specification.md see for more detail):

* A `.bed` file in Format 01.
* A `.tsv` file in Format 03.

# Config files & metadata

Config file is `config/config.Roar.yaml`. This contains paths to the `samples.csv` file, output directory, scripts, local & cluster logs and envs (conda `.yaml`'s).
Samples file present in `config/samples.csv`, and is referred to in `config/config.Roar.yaml`.
Conda environment, `workflow/envs/Roar.yaml`.



# Citation
Please use the following information to cite...

# Contact the Author
Naima Ahmed Fahmi: fnaima@knights.ucf.edu
Wei Zhang: wzhang.cs@ucf.edu
Jeongsik Yong: jyong@umn.edu
Euan McDonnell: bs14e3m@leeds.ac.uk
