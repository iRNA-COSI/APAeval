/*
 * RUN LABRAT
 */

include { LABRAT_MAKETFFASTA   } from '../modules/labrat_maketffasta' addParams( options: [:] )
include { LABRAT_RUNSALMON     } from '../modules/labrat_runsalmon'   addParams( options: [:] )

workflow RUN_LABRAT {
    take:
    ch_sample
    
    main:
    ch_sample
       .map { it -> [ it[5], it[6] ] }
       .unique()
       .set { ch_gff_fasta }

    /*
     * Run LABRAT PREPARE TF FASTA
     */
     ch_tffasta = ''
     LABRAT_MAKETFFASTA ( ch_gff_fasta )

     ch_sample
        .collect { it -> it[0] }
        .map { it.join(',') }
        .set {ch_all_samples}

     ch_sample
        .collect { it -> it[1] }
        .map { it.join(',') }
        .set { ch_all_fastq1s }

     ch_sample
        .collect { it -> it[2] }
        .map { it.unique() }
        .map { it.join(',') }
        .set { ch_all_fastq2s }

    /*
     * Run LABRAT RUNSALMON
     */
    LABRAT_RUNSALMON ( LABRAT_MAKETFFASTA.out.ch_tffasta,
                       ch_all_samples,
                       ch_all_fastq1s,
                       ch_all_fastq2s )

}

