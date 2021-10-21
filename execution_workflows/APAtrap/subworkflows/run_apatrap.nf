/*
    Run APAtrap
*/

def modules = params.modules.clone()
def run_differential = modules['final_output'].run_differential
def preprocessing    = modules['preprocessing']

include { CHECK_INPUT_PARAMS      } from '../modules/check_input_params' addParams( options: [:] )
include { CONVERT_TO_BEDGRAPH   } from '../modules/convert_to_bedgraph' addParams( options: [:] )
include { PREPROCESSING         } from '../modules/preprocessing' addParams( options: [:] )
include { IDENTIFY_DISTAL_3UTR  } from '../modules/identify_distal_3utr' addParams( options: [:] )
include { PREDICT_APA           } from '../modules/predict_apa' addParams( options: [:] )
include { CONVERT_TO_BED        } from '../modules/convert_to_bed' addParams( options: [:] )
include { DE_APA                } from '../modules/de_apa' addParams( options: [:] )
include { CONVERT_TO_TSV        } from '../modules/convert_to_tsv' addParams( options: [:] )

workflow RUN_APATRAP {
    take:
    ch_sample

    main:
    /*
        Check input params
    */
    CHECK_INPUT_PARAMS()

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
        .fromPath("$PWD/${preprocessing.genome_file}")
        .set{ ch_preprocessing_input }

    PREPROCESSING( ch_preprocessing_input )

    // rename the genome file channel for identifyDistal3UTR
    PREPROCESSING.out.ch_genome_file
        .combine(ch_3utr_input)
        .set{ ch_3utr_input }

    // If we are running differential, only the first inputs in the channels are required.
    // The input files used are in ch_sample_bedgraph_files_dir
    if ( run_differential ) {
        ch_sample_bedgraph_files_dir
            .first()
            .set{ ch_sample_bedgraph_files_dir }
    }

    /*
        Perform first step of APAtrap: identify distal 3'utr
    */
    IDENTIFY_DISTAL_3UTR( ch_sample_bedgraph_files_dir,
                          ch_3utr_input )

    // combine the bedgraph read with the inputfile for predictAPA
    ch_3utr_input
        .map{ it -> [ it[1], it[2] ] }
        .combine(IDENTIFY_DISTAL_3UTR.out.ch_predictapa_input, by:0)
        .unique()
        .set { ch_predict_apa_input }

    if ( run_differential ) {
        ch_predict_apa_input
            .first()
            .set{ ch_predict_apa_input }
    }

    /*
        Perform second step of APAtrap: infer all potential APA sites
    */
    PREDICT_APA( ch_sample_bedgraph_files_dir,
                 ch_predict_apa_input )

    // Convert PredictAPA output into identification challenges files
    if (!run_differential) {
        CONVERT_TO_BED(PREDICT_APA.out.ch_de_apa_input)
    }

    if (run_differential) {
        /*
            Perform third step of APAtrap: detect genes having significant changes
            in APA site usage between conditions
        */
        DE_APA( PREDICT_APA.out.ch_de_apa_input )

        /*
            Convert deAPA output into differential challenge tsv file
        */
        CONVERT_TO_TSV( DE_APA.out.ch_de_apa_output )
    }
}

