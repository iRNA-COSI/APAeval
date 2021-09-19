#!/usr/bin/env nextflow

// write help message here
if (params.help) {
	
	    log.info"""
	    ==============================================
	    APAeval Aptardi execution workflow
	    ==============================================
	    Usage:
	    Run the pipeline with default parameters:
	    nextflow run main.nf -profile docker --input samplesheet.csv
	    Mandatory arguments:
                --input		samplesheet
	    Flags:
                --help			Display this message
	    """.stripIndent()

	exit 1
}

// define parameters (add more parameters if necessary)
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }
ch_output_bed = params.output_bed
if (params.aptardi_model) { ch_aptardi_model = file(params.aptardi_model, checkIfExists: true) } // if not specify, the default one will be used
if (params.aptardi_scale) { ch_aptardi_scale = file(params.aptardi_scale, checkIfExists: true) } // if not specify, the default one will be used

/*
 * Check input files indicated on the samplesheet to make sure that the input files are valid input for the downstream processes
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet from ch_input
    
    output:
    path "*reformat.csv" into ch_samplesheet_reformat
    
    script:
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}

def get_sample_info(LinkedHashMap sample) {
    return [ sample.sample, sample.bam, sample.gtf, sample.fasta ]
}

// create input channel for aptardi
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { get_sample_info(it) }
    .set { ch_input }

/*
 * Run APTARDI
 */
process APTARDI {
    publishDir "${params.outdir}/aptardi", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(bam), path(gtf), path(fasta) from ch_input
    path aptardi_model from ch_aptardi_model
    path aptardi_scale from ch_aptardi_scale

    output:
    tuple val(sample), path("$sample/$output_gtf") into ch_aptardi_output

    script:
    output_gtf = "$sample"+".gtf"
    """
    aptardi --o $sample --f $fasta --r $gtf --b $bam -n model.hdf5 -t scale.pk -g $output_gtf
    """
}

/*
 * Run make output bed
 */
process MAKE_IDENTIFICATION_BED {
    tag "$sample"
    publishDir "${params.outdir}/aptardi/$sample", mode: params.publish_dir_mode //each sample has it specific output directory for the process

    input:
    val output_bed from ch_output_bed
    tuple val(sample), path(output_gtf) from ch_aptardi_output 

    output:
    path "*" into ch_bed_output

    script:
    """
    make_identification_bed.py $output_gtf $output_bed
    """
}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
