/*
 * Run CSI-UTR
 */

def modules = params.modules.clone()
def files   = modules['files']
def run_mode = modules['run_mode']
def run_identification = run_mode.run_identification
def run_differential = run_mode.run_differential

include { CREATE_INPUT_FILE } from '../modules/create_input_file'
include { CREATE_SAMPLE_INFORMATION_FILE } from '../modules/create_sample_information_file' addParams( options: [:] )
include { CSI_UTR_MAIN } from '../modules/csi_utr_main' addParams( options: [:] )
include { POSTPROCESS_IDENTIFICATION } from '../modules/postprocess_identification' addParams( options: [:] )
include { POSTPROCESS_DIFFERENTIAL } from '../modules/postprocess_differential' addParams( options: [:] )
include { PUBLISH_IDENTIFICATION_OUTPUT } from '../modules/publish_identification_output' addParams( options: [:] )
include { PUBLISH_DIFFERENTIAL_OUTPUT } from '../modules/publish_differential_output' addParams( options: [:] )

workflow RUN_CSI_UTR {
    take:
    ch_input
    ch_sample

    main:

    /*
        Move bam and bai files to an input directory
    */
    CREATE_INPUT_FILE ( ch_sample )

    /*
        Create sample information file of sample and condition names
    */
    CREATE_SAMPLE_INFORMATION_FILE ( ch_input )

    /*
        Run CSI-UTR
    */
    ch_sample
        .map { it -> [ it[0], it[3] ] }
        .set{ ch_postprocess }

    Channel
        .fromPath("${files.CSI_bed_file}")
        .set{ch_csi_bed_file}

    Channel
        .fromPath("${files.CSI_annotation_file}")
        .set{ch_csi_annotation_file}

    ch_csi_bed_file
        .combine(ch_csi_annotation_file)
        .set{ch_bed_files}
    
    CSI_UTR_MAIN (
        ch_postprocess,
        CREATE_SAMPLE_INFORMATION_FILE.out.ch_sample_information_file,
        ch_bed_files
    )

    /*
        Get identification challenge output file
    */
    if(run_identification) {
        ch_postprocess
	    .combine(CSI_UTR_MAIN.out.ch_csi_utr_outdir)
            .set { ch_postprocess_identification_input }
        POSTPROCESS_IDENTIFICATION (
            ch_postprocess_identification_input
        )
        PUBLISH_IDENTIFICATION_OUTPUT (
            POSTPROCESS_IDENTIFICATION.out.ch_identification_out
        )
    }

    /*
        Get differential challenge output file
    */
    if(run_differential) {
        POSTPROCESS_DIFFERENTIAL (
            CSI_UTR_MAIN.out.ch_csi_utr_outdir
	)
        PUBLISH_DIFFERENTIAL_OUTPUT (
            POSTPROCESS_DIFFERENTIAL.out.ch_differential_out
        )
    }
}

