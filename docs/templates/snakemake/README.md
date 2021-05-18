# {Snakemake pipeline template}

## How to create your snakemake pipeline for APAeval
Read this section, but do NOT include it your final README.
> * Copy the whole directory `snakemake` with its' subfolders to suitable location.
> * Rename `snakemake` to `[METHOD]_snakemake`, [METHOD] being the name of the tool you're building a workflow for.
> * Adjust names and contents of files that appear in this [NOTATION].
> * Place any scripts or subworkflows you're going to use into the directories `workflow/scripts` and `workflow/rules`, respectively.
> * Write your workflow in `Snakefile`.
> * Here are the [snakemake docs](https://snakemake.readthedocs.io/en/stable/index.html)
> * General good practices:
>     * Test your code with [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#best-practices).
>     * One (shell) command per rule. 
>     * If samples differ in a meaningful way (e.g. single end and paired end samples), it might be better to write subworkflows within `workflow/rules`. 
>     * If applicable, each rule has it's own conda environment. 
> * There are some shell scripts which can be used to start a snakemake run. Adjust the name of the `config file` in these scripts:
>     * `dryrun.sh`
>     * `rulegraph.sh` (Also adjust name of output `.png`)
>     * `run_local.sh`
> * Adjust this `README.md`: Delete this description section and populate the sections below with your awesome experience ;)

# [METHOD]

{Description of method, with link to publication and source code.}

## Rulegraph

{The rulegraph gives an overview of the steps of the workflow. To obtain it, adapt and run the `rulegraph.sh` script!}

![rulegraph](rulegraph.[METHOD].png)

## Input & pre-processing

{Describe input files and how they will be processed in order for the method to work. Describe how sample tables have to look like, and any other input that is needed (e.g. genome).}

## Params

{Describe parameters needed by your METHOD.}

## Output & post-processing

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
