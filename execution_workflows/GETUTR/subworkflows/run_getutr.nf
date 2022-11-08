/*
 * RUN GETUTR
 */

include { GETUTR_PROCESS  } from '../modules/getutr_process' addParams( options: [:] )
include { POSTPROCESSING_IDENTIFICATION } from '../modules/postprocessing_identification' addParams( options: [:] )
include { POSTPROCESSING_RELATIVE_QUANTIFICATION } from '../modules/postprocessing_relative_quantification' addParams( options: [:] )

workflow RUN_GETUTR {
    take:
    ch_sample
    ch_gtf
    
    main:
    ch_sample
       .map { it -> [ it[0], it[1], ch_gtf ] }
       .set { sample_bam }

    GETUTR_PROCESS(sample_bam)

    GETUTR_PROCESS.out.ch_getutr_output
        .set { ch_postprocessing_input }

    POSTPROCESSING_IDENTIFICATION(ch_postprocessing_input)

    POSTPROCESSING_RELATIVE_QUANTIFICATION(ch_postprocessing_input)
}

