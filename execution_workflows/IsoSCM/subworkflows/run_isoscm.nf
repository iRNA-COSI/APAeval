/*
 * Run IsoSCM
 */

def modules = params.modules.clone()
def run_identification = modules['run_mode'].run_identification
def run_differential = modules['run_mode'].run_differential

include { ISOSCM_ASSEMBLE } from '../modules/isoscm_assemble' addParams( options: [:] )
include { POSTPROCESS_IDENTIFICATION } from '../modules/postprocess_identification' addParams( options: [:] )

workflow RUN_ISOSCM {
     take:
     ch_star_genome_index
     ch_aligned_bam_files_dir

     main:
     ch_aligned_bam_files_dir
         .first()
         .set{ ch_aligned_bam_files_dir }

    /*
        Run IsoSCM assemble step
    */
     ISOSCM_ASSEMBLE( ch_aligned_bam_files_dir )

     if(run_identification) {
        /*
            Get identification output file
        */
         POSTPROCESS_IDENTIFICATION(ISOSCM_ASSEMBLE.out.ch_isoscm_assemble_out)
     }

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
}

