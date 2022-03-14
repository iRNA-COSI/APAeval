/*
 * RUN LABRAT
 */

include { LABRAT_MAKETFFASTA   } from '../modules/labrat_maketffasta' addParams( options: [:] )
include { SALMON_INDEX     } from '../modules/salmon_index'   addParams( options: [:] )
include { SALMON_QUANT     } from '../modules/salmon_quant'   addParams( options: [:] )
include { MAKE_QUANT_BED     } from '../modules/make_quant_bed'   addParams( options: [:] )
include { LABRAT_CALCULATEPSI   } from '../modules/labrat_calculatepsi' addParams( options: [:] )
include { MAKE_DIFFERENTIAL_TSV   } from '../modules/make_differential_tsv' addParams( options: [:] )

workflow RUN_LABRAT {
    take:
    ch_sample
    ch_sampconds
    ch_labrat_quant_dir
    conditionA
    conditionB
    run_quantification
    run_differential
    
    main:
    ch_sample
       .map { it -> [ it[3], it[4] ] }
       .unique()
       .set { ch_gff_fasta }

    /*
     * Run LABRAT PREPARE TF FASTA
     */
    ch_tffasta = ''
    LABRAT_MAKETFFASTA ( ch_gff_fasta )

    /*
     * Replace LABRAT RUNSALMON with SALMON
     */
    if (run_quantification){
        SALMON_INDEX ( LABRAT_MAKETFFASTA.out.ch_tffasta )
        ch_txfasta_idx = SALMON_INDEX.out.ch_txfasta_idx
        ch_txfasta_idx
           .combine ( ch_sample )
           .map { it -> [ it[0], it[1], it[2], it[3] ] }
           .set { ch_input_salmon_quant }
        SALMON_QUANT ( ch_input_salmon_quant )

        ch_gff_fasta
            .map { it -> it[0] }
            .set { ch_gff }

        ch_gff
            .combine( SALMON_QUANT.out.ch_salmon_quant_outputs )
            .set { ch_make_quant_bed_in }
        MAKE_QUANT_BED ( ch_make_quant_bed_in )
    }

    if (run_differential){
        if (run_quantification){
            MAKE_QUANT_BED.out.bed
                .collect()
                .set { ch_quant_done }
            SALMON_QUANT.out.ch_salmon_dir
                .map { it -> it[1] }
                .unique()
                .set { ch_salmon_dir }
            LABRAT_CALCULATEPSI ( ch_salmon_dir,
                                  ch_sampconds,
                                  conditionA,
                                  conditionB,
                                  ch_gff,
                                  LABRAT_MAKETFFASTA.out.ch_gffutilsdb,
                                  ch_quant_done )
        } else {
            ch_sample
                .map { it -> it[3] }
                .set { ch_gff }
            LABRAT_CALCULATEPSI ( ch_labrat_quant_dir,
                                  ch_sampconds,
                                  conditionA,
                                  conditionB,
                                  ch_gff,
                                  LABRAT_MAKETFFASTA.out.ch_gffutilsdb,
                                  '' )
        }
        MAKE_DIFFERENTIAL_TSV ( LABRAT_CALCULATEPSI.out.ch_pval_results,
                                conditionA, conditionB )
    }
}

