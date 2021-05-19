#!/usr/bin/env nextflow

// write help message here
if (params.help) {
	
	    log.info"""
	    ==============================================
	    APAeval pilot benmarking pipeline
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

/*
 * Check input samplesheet
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
    return [ sample.sample, sample.fastq1, sample.fastq2, sample.bam, sample.bai, sample.gff, sample.fasta, sample.bed, sample.mart_export ]
}

// create input channels for the following process(es) 
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { get_sample_info(it) }
    .into {
        ch_input_for_process1
        ch_input_for_process2
    }

if (!params.skip_process1){

    
    /*
     * Run process 1.1
     */
    process PROCESS1_1 {
        publishDir "${params.outdir}/process1", mode: params.publish_dir_mode

        input:
        path input from ch_input_for_process1

        output:
        tuple path(input_file1), val(sample_name), path(input_file2) into ch_output_process1_1

        script:
        output_dir        = "outputs"
        """
        [software] $input $output_dir
        """
    }

    // there can be multiple parts within one process
    /*
     * Run Run process 1.2
     */
    process PROCESS1_2 {
        tag "$sample"
        publishDir "${params.outdir}/process1/$sample", mode: params.publish_dir_mode //each sample has it specific output directory for the process

        input:
        tuple path(input_file1), val(sample_name), path(input_file2) from ch_output_process1_1

        output:
        path "*" into ch_output_process1_2

        script:
        output_dir = "outputs"
        """
        [software] $input_file1 $input_file2 $output_dir
        """
    }
}

// if there is a process2 then unmute the following
//if (!params.skip_process2){
//    /*
//     * Run process 2.1
//     */
//    process PROCESS2_1 {
//        publishDir "${params.outdir}/process2", mode: params.publish_dir_mode
//
//        input:
//        tuple path(input_file1), val(sample_name), path(input_file2) from ch_input_for_process2
//
//        output:
//        path "*" into ch_process2_outputs
//
//        script:
//        output_dir = "outputs"
//        """
//        [software] $input_file1 $input_file2 $output_dir
//        """
//    }
//}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
