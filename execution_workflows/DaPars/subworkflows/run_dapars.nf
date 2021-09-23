/*
 * Run DaPars
 */

def modules = params.modules.clone()
def files    = modules['files']
def mode = modules['mode']

include { PREPROCESSING         } from '../modules/preprocessing' addParams( options: [:] )
include { DAPARS_EXTRACT_3UTR   } from '../modules/dapars_extract_3utr' addParams( options: [:] )
include { CONVERT_TO_BEDGRAPH   } from '../modules/convert_to_bedgraph' addParams( options: [:] )
include { CREATE_CONFIG_FILE    } from '../modules/create_config_file' addParams( options: [:] )
include { DAPARS_MAIN           } from '../modules/dapars_main' addParams( options: [:] )
include { POSTPROCESSING        } from '../modules/postprocessing' addParams( options: [:] )

workflow RUN_DAPARS {
    take:
    ch_sample
    
    main:
    /*
        Convert gtf genome file to bed12
    */
    Channel
        .fromPath("$PWD/${files.genome_file}")
        .set{ ch_preprocessing_input }
    PREPROCESSING( ch_preprocessing_input )

    PREPROCESSING.out.ch_genome_file
        .set { ch_extract_3utr_input }

    /*
     * Run step 1 of DaPars: Generate region annotation, extract 3utr
     */
     DAPARS_EXTRACT_3UTR ( ch_extract_3utr_input )

    ch_sample
        .map { it -> [ it[0], it[1], it[2], it[3] ] }
        .set { ch_convert_to_bedgraph_input }

    /*
     * Convert sample bam files to bedgraph
     */
    CONVERT_TO_BEDGRAPH ( ch_convert_to_bedgraph_input )

    /*
     * Create config file to be used as input for step 2 of DaPars
     */
    CREATE_CONFIG_FILE (
        CONVERT_TO_BEDGRAPH.out.ch_convert_to_bedgraph_out.first(),
        DAPARS_EXTRACT_3UTR.out.ch_extracted_3utr_output
    )

    /*
     * Run step 2 of DaPars: identify the dynamic APA usages between two conditions.
     */
    DAPARS_MAIN ( CREATE_CONFIG_FILE.out.ch_dapars_input )

    /*
     * Convert DaPars output file to differential challenge output file
     */
    POSTPROCESSING ( DAPARS_MAIN.out.ch_dapars_output, mode )
}

