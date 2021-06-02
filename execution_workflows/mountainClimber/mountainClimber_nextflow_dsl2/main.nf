#!/usr/bin/env nextflow

/*
=================================================
   APAEVAL EXECUTION WORKFLOW: MOUNTAINCLIMBER
=================================================
*/

////////////////////////////////////////////////////
/* --               DECLARATIONS               -- */
////////////////////////////////////////////////////

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --              IMPORT MODULES              -- */
////////////////////////////////////////////////////

include {
    generate_bedgraphs;
} from './modules/bedtools'
include {
    prepare_reference;
    calculate_expression;
} from './modules/rsem'
include {
    prepare_index;
    align_genome;
    align_transcriptome;
} from './modules/star'
include {
//    diff_cluster;
//    diff_ru;
//    diff_segment_read_counts;
//    diff_test;
    get_junction_reads;
    merge_tus;
//    mountain_climber_cp;
//    mountain_climber_ru;
    mountain_climber_tu;
} from './modules/mountain_climber'
//include {x; y; z} from './modules/post_processing'

////////////////////////////////////////////////////
/* --            RUN MAIN WORKFLOW             -- */
////////////////////////////////////////////////////

workflow {

    take:
// TODO: set inputs

    main:

      // pre-processing
      prepare_index()

      // workflow following mountainClimber tutorial: https://github.com/gxiaolab/mountainClimber#tutorial-with-test-dataset
      align_genome()
      get_junction_reads()
      generate_bedgraphs()
      mountain_climber_tu()
      merge_tus()
      prepare_reference_rsem()
      align_transcriptome()
      calculate_expression_rsem()
//      mountain_climber_cp()
//      mountain_climber_ru()
//      diff_cluster()
//      diff_segment_read_counts()
//      diff_ru()
//      diff_test()

      // post-processing


}

////////////////////////////////////////////////////
/* --                 THE END                  -- */
////////////////////////////////////////////////////
