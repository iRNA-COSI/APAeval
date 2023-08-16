#!/usr/bin/env nextflow

// write help message here
if (params.help) {
	
	    log.info"""
	    ==============================================
	    APAeval pilot benchmarking pipeline
	    ==============================================
	    Usage:
	    Run the pipeline with default parameters:
	    nextflow run main.nf -profile docker --input samplesheet.csv
	    Mandatory arguments:
                --input   input.csv, see the example file
	    Flags:
                --help			Display this message
	    """.stripIndent()

	exit 1
}

// define parameters (add more parameters if necessary)
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }

/*
 * Step 1: simply check the input file, also output a temporary rds file
 */
process checkInput {
    tag "$inputFile"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path inputFile from ch_input
    
    output:
    path "*.rds" into ch_input_tmp
    
    script:
    """
    Rscript $baseDir/bin/check.input.R $inputFile tmp.bed.rds
    """
}


/*
 * Step 2: run the diffUTR package
 */
process runDiffUTR {
    publishDir params.outdir, mode: params.publish_dir_mode

    input:
    path inputFile from ch_input
    path tmpRDS from ch_input_tmp
    
    output:
    path "*result.tsv"

    script:
    """
    Rscript $baseDir/bin/run.diffUTR.R $inputFile $tmpRDS full.result.tsv
    # geneID,effectSize,q-val only
    cut -f1,4,12 full.result.tsv | sed '1d' > final.result.tsv
    """
}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
