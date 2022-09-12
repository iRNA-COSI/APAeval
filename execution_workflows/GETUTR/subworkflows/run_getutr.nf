/*
 * RUN LABRAT
 */

include { GETUTR_PROCESS  } from '../modules/getutr_process' addParams( options: [:] )

workflow RUN_GETUTR {
    take:
    ch_sample
    
    main:
    ch_sample
       .map { it -> [ it[0], it[1], it[2] ] }
       .set { sample_bam_gtf }

    GETUTR_PROCESS ( sample_bam_gtf )
}

