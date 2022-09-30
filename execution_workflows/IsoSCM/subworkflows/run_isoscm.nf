/*
 * Run IsoSCM assemble
 */

def modules = params.modules.clone()
def run_mode = modules['run_mode']
def run_identification = run_mode.run_identification
def run_relative_usage_quantification = run_mode.run_relative_usage_quantification
include { ISOSCM_ASSEMBLE } from '../modules/isoscm_assemble' addParams( options: [:] )
include { ISOSCM_COMPARE } from '../modules/isoscm_compare' addParams( options: [:] )
include { POSTPROCESS_IDENTIFICATION } from '../modules/postprocess_identification' addParams( options: [:] )
include { POSTPROCESS_RELATIVE_USAGE_QUANTIFICATION } from '../modules/postprocess_relative_usage_quantification' addParams( options: [:] )

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
       Get identification and/or relative usage quantification output(s)
     */
     if ( run_identification ) {
         POSTPROCESS_IDENTIFICATION( ISOSCM_COMPARE.out.ch_isoscm_compare_out )
     }

     if ( run_relative_usage_quantification ) {
         POSTPROCESS_RELATIVE_USAGE_QUANTIFICATION( ISOSCM_COMPARE.out.ch_isoscm_compare_out )
     }
}

