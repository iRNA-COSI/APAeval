/*
 * Run IsoSCM
 */

def modules = params.modules.clone()
def files   = modules['files']

include { PREPROCESS_GENOME   } from '../modules/preprocess_genome' addParams( options: [:] )

workflow RUN_ISOSCM {
    take:
    ch_sample
    
    main:
    /*
        Blabla
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

}

