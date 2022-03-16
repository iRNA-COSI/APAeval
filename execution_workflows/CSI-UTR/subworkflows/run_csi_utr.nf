/*
 * Run CSI-UTR
 */

def modules = params.modules.clone()
def files   = modules['files']

// include { CHECK_INPUT_PARAMS      } from '../modules/check_input_params' addParams( options: [:] )
include { CONVERT_GTF_TO_BED      } from '../modules/convert_gtf_to_bed' addParams( options: [:] )
include { CREATE_INPUT_FILE } from '../modules/create_input_file'
include { CREATE_SAMPLE_INFORMATION_FILE } from '../modules/create_sample_information_file' addParams( options: [:] )
include { CSI_UTR_MAIN } from '../modules/csi_utr_main' addParams( options: [:] )
// include { DAPARS_EXTRACT_3UTR     } from '../modules/dapars_extract_3utr' addParams( options: [:] )
// include { CONVERT_TO_BEDGRAPH     } from '../modules/convert_to_bedgraph' addParams( options: [:] )
// include { CREATE_CONFIG_FILE      } from '../modules/create_config_file' addParams( options: [:] )
// include { DAPARS_MAIN             } from '../modules/dapars_main' addParams( options: [:] )
// include { POSTPROCESSING          } from '../modules/postprocessing' addParams( options: [:] )

workflow RUN_CSI_UTR {
    take:
    ch_input
    ch_sample

    main:
//     /*
//         Check input params
//     */
//     CHECK_INPUT_PARAMS()

    /*
        Convert gtf genome file to bed12
    */
//     Channel
//         .fromPath("${files.genome_file}")
//         .set{ ch_gtf_genome_file }
//     CONVERT_GTF_TO_BED ( ch_gtf_genome_file )

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
//     CSI_UTR_MAIN (
//         CREATE_SAMPLE_INFORMATION_FILE.out.ch_sample_information_file,
//         CONVERT_GTF_TO_BED.out.ch_genome_file
//     )
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
        CREATE_SAMPLE_INFORMATION_FILE.out.ch_sample_information_file,
        ch_bed_files
    )

//     /*
//      * Run step 1 of DaPars: Generate region annotation, extract 3utr
//      */
//      CREATE_GENE_SYMBOL_FILE.out.ch_gene_symbol_file
//         .combine( PREPROCESSING.out.ch_genome_file )
//         .set { ch_extract_3utr_input }
//
//      DAPARS_EXTRACT_3UTR ( ch_extract_3utr_input )
//
//     ch_sample
//         .map { it -> [ it[0], it[1], it[2], it[3] ] }
//         .set { ch_convert_to_bedgraph_input }
//
//     /*
//      * Convert sample bam files to bedgraph
//      */
//     CONVERT_TO_BEDGRAPH ( ch_convert_to_bedgraph_input )
//
//     /*
//      * Create config file to be used as input for step 2 of DaPars
//      */
//     CREATE_CONFIG_FILE (
//         CONVERT_TO_BEDGRAPH.out.ch_convert_to_bedgraph_out.first(),
//         DAPARS_EXTRACT_3UTR.out.ch_extracted_3utr_output
//     )
//
//     /*
//      * Run step 2 of DaPars: identify the dynamic APA usages between two conditions.
//      */
//     DAPARS_MAIN ( CREATE_CONFIG_FILE.out.ch_dapars_input )
//
//     /*
//      * Convert DaPars output file to differential challenge output file
//      */
//      DAPARS_MAIN.out.ch_dapars_output
//         .set { ch_postprocessing_input }
//     POSTPROCESSING ( ch_postprocessing_input )
}

