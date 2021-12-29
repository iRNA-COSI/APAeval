/*
 * Run IsoSCM
 */

def modules = params.modules.clone()
def files   = modules['files']
def run_mode = modules['run_mode']

include { PREPROCESS_GENOME   } from '../modules/preprocess_genome' addParams( options: [:] )

workflow RUN_ISOSCM {
    take:
    ch_sample
    
    main:
    /*
        Check input parameters
    */
//      CHECK_INPUT_PARAMS()

    /*
        Run STAR alignment on genome files
    */
    Channel
        .fromPath("${files.gtf_genome_file}")
        .set{ ch_gtf_genome_file }

    Channel
        .fromPath("${files.fasta_genome_file}")
        .set{ ch_fasta_genome_file }

    PREPROCESS_GENOME (
            ch_gtf_genome_file,
            ch_fasta_genome_file
    )

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
}

