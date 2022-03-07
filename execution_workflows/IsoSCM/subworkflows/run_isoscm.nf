/*
 * Run IsoSCM
 */

// def modules = params.modules.clone()
//
// workflow RUN_ISOSCM {
//     take:
//     ch_sample
//
//     main:
    /*
        Check input parameters
    */
//      CHECK_INPUT_PARAMS()




//     if (run_mode.differential) {
            /*
        Merge bam files of the same condition
    */

//     SORT_BAM_FILES()

//     MERGE_BAM_FILES()

    /*
        Get bai files of the merged bam files.
        Bai files are located in the same directory as bam files
    */
//     GET_BAI_FILES()
//     }

    /*
        Run IsoSCM assemble step
    */
//     ISOSCM_ASSEMBLE()

//     if(run_mode.identification) {
        /*
            Get identification output file
        */
//         POSTPROCESS_IDENTIFICATION(ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out)
//     }

//     if(run_mode.quantification) {
        /*
        Get quantification output file
        */
//         POSTPROCESS_QUANTIFICATION(ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out)
//     }

//     if(run_mode.differential) {
    /*
        Run IsoSCM compare step
    */
//     ISOSCM_COMPARE(ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out)
    /*
        Get differential output file
    */
//     POSTPROCESS_DIFFERENTIAL(ISOSCM_COMPARE.out.ch_isoscm_compare_out)
//     }
// }

