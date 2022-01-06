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
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) } else { exit 1, "FASTA file not specified!" }
if (params.gtf)   { ch_gtf   = file(params.gtf, checkIfExists: true)   } else { exit 1, "GTF file not specified!" }
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
    return [ sample.sample, sample.bam ]
}

// create input channel for aptardi
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { get_sample_info(it) }
    .into { ch_input
            ch_stringtie_input }

if (params.use_stringtie2_gtf){
    /*
     * Generate GTF file from running STRINGTIE2
     */
    process STRINGTIE2 {
        tag "$sample"
        label 'process_medium'
        publishDir "${params.outdir}/aptardi", mode: params.publish_dir_mode

        input:
        tuple val(sample), path(bam) from ch_stringtie_input
        path gtf from ch_gtf
        path fasta from ch_fasta

        output:
        tuple val(sample), path("*.stringtie.gtf"), path(fasta) into ch_stringtie_gtf

        script:
        """
        stringtie -G $gtf -o ${sample}.stringtie.gtf $bam
        """
    }
    ch_input
       .join(ch_stringtie_gtf, by:0)
       .map { it -> [ it[0], it[1], it[2], it[3] ] }
       .set { ch_aptardi_input }
} else {
    ch_input
       .combine([ch_gtf])
       .combine([ch_fasta])
       .set { ch_aptardi_input }
}

/*
 * Run APTARDI
 */
process APTARDI {
    publishDir "${params.outdir}/aptardi", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(bam), path(gtf), path(fasta) from ch_aptardi_input
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
