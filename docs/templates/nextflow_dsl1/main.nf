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
                --input		samplesheet
	    Flags:
                --help			Display this message
	    """.stripIndent()

	exit 1
}

// define parameters (add more parameters if necessary)
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }

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
    return [ sample.sample, sample.fastq1, sample.fastq2, sample.bam, sample.bai, sample.gff, sample.fasta, sample.bed, sample.mart_export ]
}

// create input channels for the following process(es) 
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { get_sample_info(it) }
    .into {
        ch_input_for_process1_1
        ch_input_for_process1_2
        ch_input_for_process2_1
    }

// The following shows two processes: process1 and process2
// - Each of the processes can contain one step (show in process2) or multiple substeps (shown in process1)
// - The --skip_process[1/2] parameter to add flexibility to run specific processes indicated by the user.
//      Here is an example for running process1 only: nextflow main.nf --input samplesheet.csv --skip_process2

if (params.run_benchmarking_event1){

    // if the process requires a reference preprocessing step (ex. index_gff required by MISO)
    // you want to make sure that all samples in one samplesheet uses the same reference files (fasta, gff, mart_export.txt, and etc.)
    // then you can use the following nextflow function to select the unique reference (there should only be one)
    // this will process the reference for only once, which can save running time and computational resources
    ch_input_for_process1_1
       .map{ it -> it[5] }
       .unique()
       .set { ch_input_for_process1_1_gff_only }

    /*
     * Run process 1.1 (processing the reference only once)
     */
    process PROCESS1_1 {
        publishDir "${params.outdir}/process1", mode: params.publish_dir_mode

        input:
        path input from ch_input_for_process1_1_gff_only

        output:
        path "$output_dir" into ch_output_process1_1

        script:
        output_dir        = "outputs"
        """
        [software] $input $output_dir
        """
    }

    // then you want to combine the processed reference to all of the sample entries with the following
    ch_output_process1_1
       .combine ( ch_input_for_process1_2 )
       .map { it -> [ it[0], it[1], it[4], it[5] ] }
       .set { ch_input_for_process1_2_the_required_ones_only }

    /*
     * Run process 1.2 (the process will run on each sample (ex. if 7 samples, then 7 times))
     */
    process PROCESS1_2 {
        tag "$sample"
        publishDir "${params.outdir}/process1/$sample", mode: params.publish_dir_mode //each sample has it specific output directory for the process

        input:
        tuple path(processed_reference), val(sample_name), path(input_file1),path(input_file2) from ch_input_for_process1_2_the_required_ones_only 

        output:
        path "*" into ch_output_process1_2

        script:
        output_dir = "outputs"
        """
        [software] $input_file1 $input_file2 $processed_reference $output_dir
        """
    }
}

// for some processes one process is enough
if (params.run_benchmarking_event2){
    // map the inputs required for the processes
    ch_input_for_process2_1
       .map{ it -> [ it[0], it[1], it[2] ] }
       .unique()
       .set { ch_input_for_process2_1_the_required_ones_only }
    /*
     * Run process 2.1
     */
    process PROCESS2_1 {
        publishDir "${params.outdir}/process2", mode: params.publish_dir_mode

        input:
        tuple val(sample), path(input_file1), path(input_file2) from ch_input_for_process2_1_the_required_ones_only

        output:
        path "*" into ch_process2_outputs

        script:
        output_dir = "outputs"
        """
        [software] $input_file1 $input_file2 $output_dir
        """
    }
}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
