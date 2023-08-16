# APAeval utils

This directory provides Dockerfiles and corresponding 'recipes' for scripting/processing tasks that are common across workflows.

Contents:
- [Available utilities](#available-utilities)
  - [Convert GTF to BED12 (aka 'gene model') file](#convert-gtf-to-bed12-aka-gene-model-file)
  - [Convert .csv to .tsv (NO Docker)](#convert-csv-to-tsv-no-docker)
  - [Convert assessment json to .tsv (NO Docker)](#convert-assessment-json-to-tsv-no-docker)
  - [Filter assessment json (NO Docker)](#filter-assessment-json-no-docker)
- [Contributing instructions](#contributing-instructions)


## Available utilities

### Convert GTF to BED12 (aka 'gene model') file

**Description:** Convert a [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) annotation file to [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (also known as a 'gene model' file) format.

**Subdirectory name**: `gtf_to_bed12`

**Dockerhub link:** https://hub.docker.com/r/apaeval/gtftobed

**Dockerhub url (to use with Nextflow/Snakemake)**: `docker.io/apaeval/gtftobed:1.0`

**Usage (shell command to run the tool)**:
```
/gtfTobed.py --gtf <path to GTF file> --out_bed <name of the output BED file>
```

### Convert .csv to .tsv (NO Docker)

**Description:** Convert a `.csv` to a `.tsv` file using `pandas`. Useful for conversion BEFORE running an method workflow (see "Compatibility" below).

**Subdirectory name**: `csv2tsv`

**Compatibility:**
Currently no Docker image available, but the script can be run inside the `apaeval` conda environment. 

**Usage:**
```
python csv2tsv.py --csv samples.csv --tsv samples.tsv
```

### Convert assessment json to .tsv (NO Docker)

**Description:** Convert a summary workflow output json file to a `.tsv` file using `pandas`. json objects will be flattened, and each object is converted to a table row. Furthermore, the `metric_id` is split into 3 columns "metric", "window_size" and "site_set". If more than one json file is provided as input, resulting tables will be concatenated.

> NOTE: In the current implementation, only OEBs "assessment" type objects will be kept and APAeval *absolute quantification* specific parameters `window_size`, `site_set` and `metrics.metric_id` are expected to be present.This util will most likely have to be adapted for other events.

**Subdirectory name**: `metrics_json2tsv`

**Compatibility:**
Currently no Docker image available, but the script can be run inside the `apaeval` conda environment. 

**Usage:**
```
python metrics_json2tsv.py --file-list assessment1.json assessment2.json --output metrics.tsv
```

### Filter assessment json (NO Docker)

**Description:** Filter a summary workflow output json file, i.e. remove objects that belong to specified metrics or challenges. Parts of metric- or challenge names can be specified and all objects (participant, assessment, aggregation, manifest-like) containing those metrics or challenges are removed.

**Subdirectory name**: `filter_jsons`

**Compatibility:**
Currently no Docker image available, but the script can be run inside the `apaeval` conda environment. 

**Usage:**
```
python filter_jsons.py --file-list assessment1.json assessment2.json --out_prefix "filtered_" --b_metrics "union all_GT Spearman relative 25nt multi_matched FDR" --b_challenges "TE"
```
## Contributing instructions

Have code for a pre-processing / file wrangling task that you think may be general to other tools? To add a tool to the directory, we need you to create a pull request with the following tasks completed:

- Create a suitably named subdirectory under `utils` to house the subsequent files (e.g. gtf_to_bed12)
- Add a Dockerfile to the subdirectory to provide necessary dependencies and scripts. This should then be added to the `apaeval` team Dockerhub (if you don't have access, you can request as part of the pull request)
- Add description of task and how to run it to the README. This should follow the format below:


\## Header description of tool (placed under Available utilities)

(remove the escape slash when you copy and paste to interpret the hashes as titles)

**Description**: <short 1/2 sentence description of the task the utility completes>

**Dockerhub link**: `https://hub.docker.com/r/apaeval/<tool_name>`

**Subdirectory name**: `<subdirectory name under utils/ storing Dockerfile and necessary scripts>`

**Docker.io url (to use with Nextflow/Snakemake)**: <docker.io/apaeval/<tool_name:tag>

**Usage (shell command to run the tool)**:

```
python utility.py -i <input> -o <output>
```
