# APA-Scan
APA-Scan detects the potential Alternative Polyadenylation (APA) events in the downstream 3'-UTR of a gene among two different biological conditions.

# Installation
APA-Scan is a python tool which can be downloaded directly from github. Python (version 3.0 or higher) is required to be installed on users machine to run APA-Scan. It can work on Windows, Linux and Mac platform.

Users need to update the configuration.ini file and run the specific file APA-Scan.py to demonstrate the APA events.

$python3 APA-Scan.py

To generate plots for any selected 3'-UTR region, users need to run makePlots.py and provide input in a specific pattern: chromName:geneName:start-end

$python3 Make-plots.py

# User manual
User manual for APA-Scan is available on github https://github.com/compbiolabucf/APA-Scan/blob/master/APA-Scan%20User%20Manual.pdf

# Environments

**A Docker image still needs to be made for this**. All the needed packages and tools needed to run this are present in the conda `.yaml`. There is one package that couldn't be installed from conda, `apybiomart`, but it was installable from `pip` so is present in `workflow/envs/APA-Scan.py`

# Workflow (Snakemake)

Workflow is written in **Snakemake**, saved as a `workflow/Snakefile`. All the rules are grouped together under the "APA_SCAN" group, as most the rules are minimally computationally intensive (they are mostly just reformatting/file-editing) so are run with the most demanding rule, **run_apa_scan**. It runs the following rules:

1. **generate_configurationfile**: This takes as input the `blank_configuration.ini` file and replaces several placeholders with the desired input/output directories and files, outputting a `configuration.ini` file. Note that APA-Scan apparently needs the absolute directory path for this to work, but the rule works by getting this using the python command `os.getwd()` so that it generalises to any system.
2. **run_apa_scan**: This takes as input the `configuration.ini` file and runs `APA-Scan.py`. It outputs a `.xlsx` file.
3. **convert_to_tsv**: This uses `xlsx2csv` to convert the `.xlsx` output of `APA-Scan.py` to a `.tsv` file, with tabs as the delimiter.
4. **convert_gene_ids**: This takes as input the `.tsv` file produced from the previous rule and maps the gene names to ENSEMBL GENE ID, then reformats and saves as both a `.bed` (format 01) or a `.tsv` (format 03) as output. See https://github.com/iRNA-COSI/APAeval/blob/main/execution_workflows/execution_output_specification.md for output formats. It runs the script `workflow/scripts/convert_geneids.py` which runs the above.

It is recommended to run the following command with snakemake, to force the `configuration.ini` file to be recreated anew for each run.

```
snakemake --forcerun generate_configurationfile --snakefile workflow/Snakefile
```

# Input files

This workflow takes as input:

* At least 2 `.bam` files split between 2 different conditions, generally 1 treatment, 1 control. These should be split into 2 directories corresponding to the different conditions. So for the test data, the files would be split as `srsf3/siSrsf3_R1_100genes.bam` and `control/siControl_R1_100genes.bam`.
* A genome `.fa` file for the genome of the organism being tested. For the test data this should be `mm10.fa`, ie for the *mus musculus* genome.
* An annotation file should match the genome file and can be downloaded from UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables, see manual for more details). For the data, this is present in the directory as `mm10.ncbiRefSeq.txt`.
* A `blank_configuration.ini` file, that is processed by the workflow into the `configuration.ini` file that is input into `APA-Scan.py` (See above section).

# Output files

The workflow outputs 2 formats (https://github.com/iRNA-COSI/APAeval/blob/main/execution_workflows/execution_output_specification.md see for more detail):

* A `.bed` file in Format 01.
* A `.tsv` file in Format 03.

# Config files & metadata

Config file is `config/config.APA-Scan.yaml`. This contains paths to the `samples.csv` file, output directory, scripts, local & cluster logs and envs (conda `.yaml`'s).
Samples file present in `config/samples.csv`, and is referred to in `config/config.APA-Scan.yaml`.
Conda environment, `workflow/envs/APA-Scan.yaml`.

# Notes on APA-Scan.py

I had to manually edit the  `APA-Scan.py` and `methods.py` scripts myself as they didn't seem to work without making a few adjustments. These are commented in the files themselves and marked with "~ EM" to indicate that I (Euan McDonnell) made the edits.



# Citation
Please use the following information to cite...

# Contact the Author
Naima Ahmed Fahmi: fnaima@knights.ucf.edu
Wei Zhang: wzhang.cs@ucf.edu
Jeongsik Yong: jyong@umn.edu
Euan McDonnell: bs14e3m@leeds.ac.uk
