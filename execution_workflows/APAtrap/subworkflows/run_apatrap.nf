/*
 * Run APAtrap
 */

include { PREPROCESSING } from '../modules/preprocessing' addParams( options: [:] )

workflow RUN_APATRAP {
    take:
    ch_sample

    main:
    ch_sample
        .map { it -> [ it[3], it[4] ] }
        .unique()
        .set { ch_preprocess }

    /*
        Preprocess the input file: convert from bam to bedgraph
    */
    PREPROCESSING( ch_preprocess )

    /*
        Perform first step of APAtrap: identify distal 3'utr
    */
    ch_identify_3utr_input = PREPROCESSING.out.ch_preprocessing_output
    IDENTIFY_DISTAL_3UTR( ch_identify_3utr_input )
}

