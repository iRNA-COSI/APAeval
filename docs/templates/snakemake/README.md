# {Snakemake pipeline template}

## How to create your snakemake pipeline for APAeval
Read this section, but do NOT include it your final README.
> * Copy the whole directory `snakemake` with its' subfolders to suitable location.
> * Rename `snakemake` to `[METHOD]_snakemake`, [METHOD] being the name of the tool you're building a workflow for.
> * Adjust names and contents of files that appear in this [NOTATION].
> * Place any scripts or subworkflows you're going to use into the directories `workflow/scripts` and `workflow/rules`, respectively.
> * Write your workflow in `Snakefile`.
> * Here are the [snakemake docs](https://snakemake.readthedocs.io/en/stable/index.html)
> * The pipeline's config.yaml file should contain 'True/False' flags to run each benchmarking event ('run_identification', 'run_quantification', 'run_differential') with which the tool is compatible.
>     * Wherever possible, these flags should be compatible with one another (i.e. 'switching on' any combination of flags will produce a valid workflow).
>     * If workflows are incompatible with one another (e.g. 'run_differential' & 'run_quantification' require different workflows), this should be clearly stated in the README.
>     * If a tool does not output files compatible with a challenge, this should be clearly stated in the README. The  corresponding config.yaml flag should be set to False and enforced in the Snakefile (there is commented code in the template Snakefile to do this).
> * General good practices:
>     * Test your code with [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#best-practices).
>     * One (shell) command per rule.
>     * If samples differ in a meaningful way (e.g. single end and paired end samples), it might be better to write subworkflows within `workflow/rules`.
>     * Each rule has it's own [Docker container](https://www.docker.com/resources/what-container). All containers should be pushed to the [APAeval team Dockerhub](https://hub.docker.com/u/apaeval)
> * There are some shell scripts which can be used to start a snakemake run. Adjust the name of the `config file` in following scripts:
>     * `dryrun.sh`
>     * `rulegraph.sh` (Also adjust name of output `.png`)
>     * `run_local.sh`
>     * Both `dryrun.sh` and `rulegraph.sh` will execute successfully with the current template.
> * Check out the pilot benchmark at `tests/pilot_benchmark/snakemake` for a running example. It illustrates how the execution workflow can look in practice.
> * Adjust this `README.md`: Delete this description section and populate the sections below with your awesome experience ;)

# [METHOD] {Method name as specified in Algorithm table.}

## Benchmarking challenges

{Bullet point list of the benchmarking challenges with which this workflow is compatible and produces challenge output files}

## Rulegraph

{The rulegraph gives an overview of the steps of the workflow. To obtain it, adapt and run the `rulegraph.sh` script!}

![rulegraph](rulegraph.[METHOD].png)

## Input & pre-processing

{Describe input files and how they will be processed in order for the method to work. Describe how sample tables have to look like, and any other input that is needed (e.g. genome).

Note: A single sample table format should be sufficient for all workflows (different columns will be used based on the run mode). If not possible, clearly state that a different sample table format is required for different workflows and provide an example of each possible format.}

## Params

{Describe parameters needed by your METHOD.}

{Describe the parameters for running different modes `run_identification`, `run_quantification`, `run_differential`}

## Output & post-processing

{Describe output files and postprocessing steps if necessary.}

## Notes

{Notes about the METHOD.
e.g. Did you have to adjust the method's source code/recommended running instructions?
}
