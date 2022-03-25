/*
	Generate bam index file for input bam file
*/

include { GENERATE_BAI_FILE } from '../modules/generate_bai_file' addParams( options: [:] )

workflow GENERATE_BAM_INDEX_FILE {
	take:
	ch_sample

	main:
        ch_sample
	    .map { it -> [ it[0], it[2], it[1] ] }
            .set{ ch_generate_bai_file_input }

	GENERATE_BAI_FILE (
	    ch_generate_bai_file_input
	)

	GENERATE_BAI_FILE.out.ch_bam_files
	    .set{ ch_bam_files }

	emit:
	ch_bam_files
}
        
