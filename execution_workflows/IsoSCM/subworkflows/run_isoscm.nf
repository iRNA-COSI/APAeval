/*
 * Run IsoSCM assemble
 */

def modules = params.modules.clone()
def run_identification = modules['run_mode'].run_identification
def run_relative_usage_quantification = modules['run_mode'].run_relative_usage_quantification
include { ISOSCM_ASSEMBLE } from '../modules/isoscm_assemble' addParams( options: [:] )
include { POSTPROCESS_IDENTIFICATION } from '../modules/postprocess_identification' addParams( options: [:] )
include { ISOSCM_COMPARE } from '../modules/isoscm_compare' addParams( options: [:] )
include { POSTPROCESS_RELATIVE_USAGE_QUANTIFICATION } from '../modules/postprocess_relative_usage_quantification' addParams( options: [:] )

workflow RUN_ISOSCM {
     take:
     ch_aligned_bam_files_dir

     main:
     /*
        Run IsoSCM assemble step for each sample
     */
     ISOSCM_ASSEMBLE( ch_aligned_bam_files_dir )

     ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out
         .map { it -> [ it[0], it[1] ] }
         .set{ ch_postprocess_identification_in }
     
     /*
        Get identification output file
     */
     if ( run_identification ) {
         POSTPROCESS_IDENTIFICATION( ch_postprocess_identification_in ) 
     }

    /*
        Get relative usage quantification output file
    */
    if ( run_relative_usage_quantification ) {
        ISOSCM_COMPARE( ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out_xml )
        POSTPROCESS_RELATIVE_USAGE_QUANTIFICATION( ISOSCM_COMPARE.out.ch_isoscm_compare_out )
    }
}

