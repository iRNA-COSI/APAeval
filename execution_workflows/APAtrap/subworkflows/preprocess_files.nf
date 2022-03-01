/*
    Run APAtrap
*/

def modules = params.modules.clone()
def preprocessing = modules['preprocessing']

include { CONVERT_TO_BEDGRAPH   } from '../modules/convert_to_bedgraph' addParams( options: [:] )
include { PREPROCESSING         } from '../modules/preprocessing' addParams( options: [:] )

workflow PREPROCESS_FILES {
    take:
    ch_sample

    main:

    // get the bam and bai files
    ch_sample
        .map { it -> [ it[0], it[1], it[2], it[3] ] }
        .unique()
        .set { ch_convert_to_bedgraph_input }


    /*
        Preprocess the input file: convert from bam to bedgraph
    */
    CONVERT_TO_BEDGRAPH( ch_convert_to_bedgraph_input )

    CONVERT_TO_BEDGRAPH.out.ch_sample_bedgraph_files_dir
        .set{ ch_sample_bedgraph_files_dir }

    CONVERT_TO_BEDGRAPH.out.ch_3utr_input
        .set{ ch_3utr_input }

    /*
        Convert gtf genome file to bed12
    */
    Channel
        .fromPath("${preprocessing.genome_file}")
        .set{ ch_preprocessing_input }

    PREPROCESSING( ch_preprocessing_input )

    // rename the genome file channel for identifyDistal3UTR
    PREPROCESSING.out.ch_genome_file
        .combine(ch_3utr_input)
        .set{ ch_3utr_input }

    emit:
    ch_preprocessing_input
    ch_sample_bedgraph_files_dir
    ch_3utr_input
}

