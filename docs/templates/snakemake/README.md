# Snakemake pipeline templates



## How to create own pipeline

* Copy `snakemake` to suitable location.
* Rename `snakemake` to `[METHOD]_snakemake`.
* Adjust names and contents of files that appear in this [NOTATION].
  * directories `workflow/scripts` and `workflow/rules` are optional.
* Write execution rules in `Snakefile`.
  > 1. Test your code with [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#best-practices).
  > 2. One (shell) command per rule. 
  > 3. If samples differ in meaningful way (e.g. single end and paired end samples), it might be better to write subworkflows within `workflow/rules`. 
  > 4. If applicable, each rule has an own conda environment. 
  
* Adjust configuration file name in run scripts:
  * `dryrun.sh`
  * `rulegraph.sh`
  * `run_local.sh`
* Adjust this `README.md` and populate the sections below.

# [METHOD]

{Description of method, with link to publication and source code.}

## Rulegraph

{The rulegraph gives an overview of the method. Place the image below.}

![rulegraph](rulegraph.[METHOD].png)

## Input

{Describe sample table and any other input (e.g. genome).}

## Params

{Describe parameters needed by your METHOD.}

## Output

{Describe method output and postprocessing steps if necessary.}

> Naming convention for output files:
> [IDENTIFIER]_[METHOD].[OUTFORMAT]

## Notes

{Notes about the METHOD. 
e.g. Did you have to adjust the method's soure code?
}

## Verdict

* How easy is it to install?
* How extensivly documented is the method?
* If a test is supplied, does it work?
* How easy is it to run?
* Any other experiences, remarks?
