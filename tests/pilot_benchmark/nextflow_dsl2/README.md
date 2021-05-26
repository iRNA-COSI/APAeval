# APAeval nextflow pilot run

### Steps to run this:
 - To run this first download the required reference files from [here](https://drive.google.com/drive/folders/1MUMilzaqef9u0sjScxzgi0JPTKAFQq_T?usp=sharing)
 - replace "path_to" in samplesheet_example_files.csv with the path to the reference files downloaded from the Google Drive link
 - check the path to this directory with `pwd` and replace the `[pwd]` in samplesheet_example_files.csv with the path from the `pwd` command
 - Before running `QAPA`, you have to install `QAPA`, and the guide to install is [here](https://github.com/morrislab/qapa). Make sure that QAPA is on your $PATH
 - Then, you are good to run the pilot benchmark nextflow pipeline with `QAPA`
```
nextflow main.nf --input samplesheet_example_files.csv
```
