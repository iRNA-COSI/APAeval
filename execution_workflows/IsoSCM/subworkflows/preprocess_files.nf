/*
	Preprocess files for IsoSCM by adding XS tag in BAM files to be compatible with IsoSCM
*/

def modules = params.modules.clone()
def files   = modules['files']
def run_star_genome_generate = modules['run_mode'].run_star_genome_generate

include { STAR_GENOME_GENERATE } from '../modules/star_genome_generate' addParams( options: [:] )
include { BAM_TO_FASTQ         } from '../modules/bam_to_fastq' addParams( options: [:] )
include { STAR_ALIGNMENT       } from '../modules/star_alignment' addParams( options: [:] )
include { GENERATE_BAI_FILE    } from '../modules/generate_bai_file' addParams( options: [:] )

workflow PREPROCESS_FILES {
        take:
        ch_sample
        
        main:
  
        /*
            Generate genome index with STAR if run_star_genome_generate is set to true
        */
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

            STAR_GENOME_GENERATE.out.ch_star_genome_indicator
		.set{ ch_star_genome_indicator }
            
            STAR_GENOME_GENERATE.out.ch_star_index_dir
                .set{ ch_star_genome_index }
        }
        
        /*
            If run_star_genome_generate is set to false,
            obtain genome index from the path provided by the user
        */
        else {
            ch_star_genome_indicator = Channel.of("done")
        
            Channel
                .fromPath("${files.star_genome_index}")
                .set{ ch_star_genome_index }
        }
          
        /*
            Convert input BAM file to FASTQ for STAR alignment
        */
        ch_sample
            .map { it -> [ it[0], it[1], it[3] ] }
            .unique()
            .set { ch_bam_to_fastq_input }
 
        BAM_TO_FASTQ( ch_bam_to_fastq_input )
 
        /*
            Prepare STAR ailgnment inputs
        */
        ch_sample
            .map { it -> [ it[0], it[2], it[3] ] }
            .unique()
            .set { ch_star_alignment_input }
        
        ch_star_alignment_input
	    .combine(BAM_TO_FASTQ.out.ch_fastq_files, by: 0)
            .set { ch_star_alignment_input }
       
        ch_star_alignment_input
            .combine(ch_star_genome_index)
            .set { ch_star_alignment_input }

        ch_star_alignment_input
            .combine(ch_star_genome_indicator)
            .set { ch_star_alignment_input }
            
       /*
           Run STAR alignment to obtain BAM file with XS tag
       */
       STAR_ALIGNMENT (
            ch_star_alignment_input
       )	
      
       /*
           Generate index file for the BAM file outputted by STAR alignment
       */
       GENERATE_BAI_FILE(
           STAR_ALIGNMENT.out.ch_aligned_bam_files
       )

         
       GENERATE_BAI_FILE.out.ch_bam_files
           .set{ ch_aligned_bam_files }

       emit:
       ch_aligned_bam_files
}
