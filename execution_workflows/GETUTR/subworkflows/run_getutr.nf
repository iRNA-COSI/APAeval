/*
 * RUN LABRAT
 */

include { GETUTR_PROCESS  } from '../modules/getutr_process' addParams( options: [:] )

workflow RUN_GETUTR {
    take:
    ch_sample
    ch_sampconds
    ch_labrat_quant_dir
    
    main:
    ch_sample
       .map { it -> [ it[3], it[4] ] }
       .unique()
       .set { ch_gff_fasta }

    ch_tffasta = ''
    GETUTR_PROCESS ( ch_gff_fasta )

}

