# Q1: Compute resources benchmark


## Synopsis

The goal is to assess the compute performance of an individual workflow and put in perspective of input sample size.
As the implementation and workflow language can vary from execution workflow to execution workflow, the benchmark should focus on single threaded performance and report the corresponding CPU time.
Similarly, the maximum memory usage should be reported.
The total disk space usage can also be reported. 

The three performance metrics can be plotted against the sample input size in order to assess scaling effects.


## General info

* **Challenge:** Quantification
* **Identifier:** Q1

## Inputs

| # | Format | Link | Example data |
| --- | --- | --- | --- |
| 1 | json | [Specification][spec-json] | [Link][in1] |
| 2 | json | [Specification][spec-json] | [Link][in2] |
| 3 | json | [Specification][spec-json] | [Link][in3] |

### Additional info inputs


#### json 1

The Input 1 is obtained by running the execution workflow with the `time` command in GNU version.
It includes the following entries:
* `runtime_sec`: float of user CPU time in seconds.
* `max_mem_kib`: integer of maximum memory usage in Kbytes as identified with `MaxRSS`.

Example wrapping of command `my_tool` to obtain desired json Input 1 (= output of workflow):
```bash
\time -f "{\"runtime_sec\": %U, \"max_mem_kib\": %M}" -o input1.json my_tool
```

#### json 2

Input 2 reports the the disk size usage of the output directory, without snakemake or nextflow specific directories. For example, `.snakemake` or `workflow` are excluded from the disk size usage.
It includes the following entries:
* `dir_usage_mib`: integer of disk space usage in MBytes.

Example for obtaining the size of OUTDIR
```bash
du -d 1 -BM OUTDIR 
```

#### json 3

Input 3 reports the disk size usage of the sample input file, namely the FASTQ (.fq) raw file size.
The json file includes the following entries:
* `sample_size_mib`: integer of sample FASTQ (.fq) file size in MBytes.

Example for obtaining the file size of INPUTSAMPLE.fq
```bash
du -BM INPUTSAMPLE.fq
```

## Outputs

| # | Format | Link | Example data |
| --- | --- | --- | --- |
| 1 | json | [Specification][spec-json] | [Link][out1] |


### Additional info outputs

#### json 1

The output json 1 file collects all input files into one json file.

## Metrics

| # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | CPU time for sample input file size | sec | json 1 and 3 | Can be directly obtained | matrix of sample size and CPU time | Can be plotted as 2D graph (input size versus CPU time) |
| 2 | Max memory usage for sample input file size | float | json 1 and 3 | Can be directly obtained | matrix of sample size and max memory | Can be plotted as 2D graph (input size versus max memory usage) |
| 3 | Disk space usage for sample input file size | float | json 2 and 3 | Can be directly obtained | matrix of sample size and disk space usage | Can be plotted as 2D graph (input size versus disk space usage) |


[//]: # (References)

[short-hand-ref]: <https://my-url-target.edu>
[in1]: ./example_files/input1.json
[in2]: ./example_files/input2.json
[in3]: ./example_files/input3.json
[out1]: ./example_files/output1.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
