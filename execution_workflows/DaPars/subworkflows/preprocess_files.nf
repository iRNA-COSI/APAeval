/*
    Run DaPars preprocessing steps
*/

def modules = params.modules.clone()
def files   = modules['files']

include { CHECK_INPUT_PARAMS      } from '../modules/check_input_params' addParams( options: [:] )
include { PREPROCESSING           } from '../modules/preprocessing' addParams( options: [:] )
include { CREATE_GENE_SYMBOL_FILE } from '../modules/create_gene_symbol_file' addParams( options: [:] )
include { DAPARS_EXTRACT_3UTR     } from '../modules/dapars_extract_3utr' addParams( options: [:] )
include { CONVERT_TO_BEDGRAPH     } from '../modules/convert_to_bedgraph' addParams( options: [:] )

workflow PREPROCESS_FILES {
    take:
    ch_sample

    main:
    /*
        Check input params
    */
    CHECK_INPUT_PARAMS()

    /*
        Convert gtf genome file to bed12
    */
    Channel
        .fromPath("${files.genome_file}")
        .set{ ch_gtf_genome_file }
    PREPROCESSING ( ch_gtf_genome_file )

    /*
        Create gene symbol file from gtf genome file
    */
    CREATE_GENE_SYMBOL_FILE ( ch_gtf_genome_file )

    /*
     * Run step 1 of DaPars: Generate region annotation, extract 3utr
     */
     CREATE_GENE_SYMBOL_FILE.out.ch_gene_symbol_file
        .combine( PREPROCESSING.out.ch_genome_file )
        .set { ch_extract_3utr_input }

     DAPARS_EXTRACT_3UTR ( ch_extract_3utr_input )

    ch_sample
        .map { it -> [ it[0], it[1], it[2], it[3] ] }
        .set { ch_convert_to_bedgraph_input }

    /*
     * Convert sample bam files to bedgraph
     */
    CONVERT_TO_BEDGRAPH ( ch_convert_to_bedgraph_input )

    /*
    * Rename output files
    */
    CONVERT_TO_BEDGRAPH.out.ch_convert_to_bedgraph_out
        .set{ ch_convert_to_bedgraph_out }

    DAPARS_EXTRACT_3UTR.out.ch_extracted_3utr_output
        .set{ ch_extracted_3utr_out }

    emit:
    ch_convert_to_bedgraph_out
    ch_extracted_3utr_out

}