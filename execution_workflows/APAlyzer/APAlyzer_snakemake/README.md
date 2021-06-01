# {Snakemake pipeline template}

## How to create your snakemake pipeline for APAeval
Read this section, but do NOT include it your final README.
> * Copy the whole directory `snakemake` with its' subfolders to suitable location.
> * Rename `snakemake` to `[METHOD]_snakemake`, [METHOD] being the name of the tool you're building a workflow for.
> * Adjust names and contents of files that appear in this [NOTATION].
> * Place any scripts or subworkflows you're going to use into the directories `workflow/scripts` and `workflow/rules`, respectively.
> * Write your workflow in `Snakefile`.
> * Here are the [snakemake docs](https://snakemake.readthedocs.io/en/stable/index.html)
> * Create the [conda](https://docs.conda.io/en/latest/) environment named `snakemake` with `conda env create -f snakemake.yaml`. Alternatively, create an own virtual environment, but ensure to use same versions as in `snakemake.yaml`.
> * General good practices:
>     * Test your code with [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#best-practices).
>     * One (shell) command per rule. 
>     * If samples differ in a meaningful way (e.g. single end and paired end samples), it might be better to write subworkflows within `workflow/rules`. 
>     * If applicable, each rule has it's own [Docker container](https://www.docker.com/resources/what-container) or conda environment.
> * There are some shell scripts which can be used to start a snakemake run. Adjust the name of the `config file` in these scripts:
>     * `dryrun.sh`
>     * `rulegraph.sh` (Also adjust name of output `.png`)
>     * `run_local.sh`
> * Check out the pilot benchmark at `tests/pilot_benchmark/snakemake` for a running example. It illustrates how the execution workflow can look in practice. 
> * Adjust this `README.md`: Delete this description section and populate the sections below with your awesome experience ;)


# APAlyzer
APAlyzer utilizes the PAS (polyadenylation sites) collection in the PolyA_DB database
(http://polya-db.org/polya_db/v3/) to examine APA (alternative polyadenylation) events in all genic regions,
including 3′UTRs and introns.

Publication: https://academic.oup.com/bioinformatics/article/36/12/3907/5823886
Source code: https://bioconductor.org/packages/release/bioc/html/APAlyzer.html
Software manual: https://bioconductor.org/packages/release/bioc/manuals/APAlyzer/man/APAlyzer.pdf
Github repo: https://github.com/RJWANGbioinfo/APAlyzer
Github analysis example: https://github.com/RJWANGbioinfo/APAlyzer#complete-analysis-example-apa-analysis-in-mouse-testis-versus-heart

## Input & pre-processing
* Bam files
* PAS reference regions

* bamifiles
download_testbam()
flsall <- dir(getwd(),".bam")
flsall<-paste0(getwd(),'/',flsall)
names(flsall)<-gsub('.bam','',dir(getwd(),".bam"))


* PAS reference
library(repmis)
URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
file="mm9_REF.RData"
source_data(paste0(URL,file,"?raw=True"))
PASREF=REF4PAS(refUTRraw,dfIPA,dfLE)
UTRdbraw=PASREF$UTRdbraw
dfIPA=PASREF$dfIPA
dfLE=PASREF$dfLE

* RELATIVE UTR AND INTRON EXPRESSION
UTR_APA_OUT=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="invert")
IPA_OUT=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="invert", nts=1)

* DIFF USAGE

############# 3utr APA #################
sampleTable = data.frame(samplename = 
c('Heart_rep1',
'Heart_rep2',
'Heart_rep3',
'Heart_rep4',
'Testis_rep1',
'Testis_rep2',
'Testis_rep3',
'Testis_rep4'),
condition = c(rep("Heart",4),
rep("Testis",4)))
						
test_3UTRAPA=APAdiff(sampleTable,UTR_APA_OUT, 
                        conKET='Heart',
                        trtKEY='Testis',
                        PAS='3UTR',
                        CUTreads=5)
						
table(test_3UTRAPA$APAreg)						
APAVolcano(test_3UTRAPA, PAS='3UTR', Pcol = "pvalue", plot_title='3UTR APA')
APABox(test_3UTRAPA, xlab = "APAreg", ylab = "RED", plot_title = NULL)

############# IPA #################
test_IPA=APAdiff(sampleTable,
                        IPA_OUT, 
                        conKET='Heart',
                        trtKEY='Testis',
                        PAS='IPA',
                        CUTreads=5)
table(test_IPA$APAreg)	
APAVolcano(test_IPA, PAS='IPA', Pcol = "pvalue", plot_title='IPA')
APABox(test_IPA, xlab = "APAreg", ylab = "RED", plot_title = NULL)

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

## Rulegraph

{The rulegraph gives an overview of the steps of the workflow. To obtain it, adapt and run the `rulegraph.sh` script!}

![rulegraph](rulegraph.[METHOD].png)