 {nextflow pipeline template}

## How to create your nextflow pipeline for APAeval
Read this section, but do NOT include it your final README.
> * Copy the whole directory `nextflow` with its' subfolders to suitable location.
> * Rename `nextflow` to `[METHOD]_nextflow`, [METHOD] being the name of the tool you're building a workflow for.
> * Adjust names and contents of files that appear in this [NOTATION].
> * Place any scripts or subworkflows you're going to use into the directories `workflow/scripts` and `workflow/rules`, respectively.
> * Write your workflow in `main.nf`.
> * Here are the [nextflow docs](https://www.nextflow.io/docs/latest/getstarted.html)
> * General good practices:
>     * Test your code with [nextflow main.nf --input <input_file> --other_parameters].
>     * Initiate your parameters in `nextflow.config`.
>     * Create a samplesheet.csv file so that all input files and reference files are organized and can be checked whether each file inputted is valid for running downstream process(es).
>     * If samples differ in a meaningful way (e.g. single end and paired end samples), it will be great to adjust part of the command for running a software (see example [here](https://github.com/yuukiiwa/APAeval/blob/856b8a9455f8352037e5b9a202e52b0c980d90f8/tests/pilot_benchmark/nextflow/main.nf#L156)). 
>     * If applicable, software has it's own docker container.
>     * Adjust this `README.md`: Delete this description section and populate the sections below with your awesome experience ;)

# [METHOD]

{Description of method, with link to publication and source code.}

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

