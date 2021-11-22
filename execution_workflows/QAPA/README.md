# APAeval nextflow pilot run

### Steps to run this:
 - To run this first download the required reference files from [here](https://drive.google.com/drive/folders/1MUMilzaqef9u0sjScxzgi0JPTKAFQq_T?usp=sharing)
 - replace "path_to" in samplesheet_example_files.csv with the path to the reference files downloaded from the Google Drive link
 - check the path to this directory with `pwd` and replace the `[pwd]` in samplesheet_example_files.csv with the path from the `pwd` command
 - Before running `QAPA`, you have to install `QAPA`, and the guide to install is [here](https://github.com/morrislab/qapa). Make sure that QAPA is on your $PATH
 - Then, you are good to run the pilot benchmark nextflow pipeline with `QAPA`

### Docker
This workflow uses docker containers. To run with docker, make sure that docker is installed and running (e.g. to ensure docker is running, run the command docker --help and a help message should be printed). If running with Docker, please indicate `-profile docker`.
The docker container for QAPA is available here: https://hub.docker.com/r/apaeval/qapa

### Singularity
To run with singularity, please indicate `-profile singularity`.

### To run the pipeline
```
nextflow main.nf --input samplesheet_example_files.csv -profile [docker/singularity]
```
