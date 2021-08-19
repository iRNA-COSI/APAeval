/*
 * Run DaPars
 */

include { DAPARS_EXTRACT_3UTR   } from '../modules/dapars_extract_3utr' addParams( options: [:] )

workflow RUN_DAPARS {
    take:
    ch_sample
    
    main:
    ch_sample
       .map { it -> [ it[4], it[5] ] }
       .unique()
       .set { ch_extract_3utr_input }

    /*
     * Run step 1 of DaPars: Generate region annotation, extract 3utr
     */
     DAPARS_EXTRACT_3UTR ( ch_extract_3utr_input )

    ch_sample
        .map { it -> [ it[2], it[3] ] }
        .unique()
        .set { ch_convert_to_bedgraph_input }
    /*
     * Convert sample bam files to bedgraph
     */
     CONVERT_TO_BEDGRAPH ( ch_convert_to_bedgraph_input )

    CHECK_INPUTS()

    /*
     * Create config file to be used as input for step 2 of DaPars
     */
    //CREATE_CONFIG_FILE ()
}

