/*
 * Run IsoSCM assemble
 */

def modules = params.modules.clone()

include { ISOSCM_ASSEMBLE } from '../modules/isoscm_assemble' addParams( options: [:] )
include { ISOSCM_COMPARE } from '../modules/isoscm_compare' addParams( options: [:] )
include { POSTPROCESS_IDENTIFICATION } from '../modules/postprocess_identification' addParams( options: [:] )

workflow RUN_ISOSCM {
     take:
     ch_aligned_bam_files_dir

     main:
     /*
        Run IsoSCM assemble step for each sample
     */
     ISOSCM_ASSEMBLE( ch_aligned_bam_files_dir )

     /*
        Run IsoSCM compare step for each sample
     */
     ISOSCM_COMPARE( ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out_xml )

     /*
        Get identification output file
     */
     ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out
         .map { it -> [ it[0], it[1] ] }
         .set{ ch_postprocess_identification_in }
     POSTPROCESS_IDENTIFICATION( ch_postprocess_identification_in ) 
}

