/*
	Preprocess files for IsoSCM
*/

def modules = params.modules.clone()
def files   = modules['files']
def run_star_genome_generate = modules['run_mode'].run_star_genome_generate

include { STAR_GENOME_GENERATE } from '../modules/star_genome_generate' addParams( options: [:] )
include { STAR_ALIGNMENT       } from '../modules/star_alignment' addParams( options: [:] )

workflow PREPROCESS_FILES {
        take:
        ch_sample
        
        main: 
        if( run_star_genome_generate ) {
            Channel
                .fromPath("${files.gtf_genome_file}")
                .set{ ch_gtf_genome_file }

            Channel
                .fromPath("${files.fasta_genome_file}")
                .set{ ch_fasta_genome_file }

            /*
                Run STAR alignment on genome files
            */
            STAR_GENOME_GENERATE (
                ch_gtf_genome_file,
                ch_fasta_genome_file
            )
            ch_star_genome_indicator = STAR_GENOME_GENERATE.out.ch_star_genome_index
        }
        
        if ( run_star_genome_generate == false ) {
            ch_star_genome_indicator = Channel.of("done")
        }
     
        ch_sample
            .map { it -> [ it[0], it[1], it[2], it[3] ] }
            .unique()
            .set { ch_sample }
        
        Channel
            .fromPath("${files.star_genome_index}")
            .set{ ch_star_genome_index }
        
        ch_sample
	    .combine(ch_star_genome_index)
            .set{ ch_star_alignment_input }

        ch_star_alignment_input
            .combine( ch_star_genome_indicator )
            .set{ ch_star_alignment_input }
        
        STAR_ALIGNMENT (
            ch_star_alignment_input
        )	
      
       STAR_ALIGNMENT.out.ch_aligned_bam_files_dir
           .set{ ch_aligned_bam_files_dir }
       
       emit:
       ch_star_genome_index
       ch_aligned_bam_files_dir
}
