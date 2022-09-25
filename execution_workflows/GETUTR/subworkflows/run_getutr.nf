/*
 * RUN GETUTR
 */

include { GETUTR_PROCESS  } from '../modules/getutr_process' addParams( options: [:] )

workflow RUN_GETUTR {
    take:
    ch_sample
    ch_gtf
    
    main:
    ch_sample
       .map { it -> [ it[0], it[1], ch_gtf ] }
       .set { sample_bam }

    GETUTR_PROCESS ( sample_bam )
}

