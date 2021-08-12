/*
    Run APAtrap
*/

def workflow_option = params.workflow.clone()
def run_differential = workflow_option['run_differential']

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
    // get the bam and bai files
    ch_sample
        .map { it -> [ it[0], it[1], it[2] ] }
        .unique()
        .set { ch_convert_to_bedgraph_input }


    /*
        Preprocess the input file: convert from bam to bedgraph
    */
    CONVERT_TO_BEDGRAPH( ch_convert_to_bedgraph_input )

    CONVERT_TO_BEDGRAPH.out.ch_bedgraph
        .set { ch_preprocessing_input }

    /*
        Ensure that file has leading chr
    */
    PREPROCESSING( ch_preprocessing_input )

    // combine the genome file with the input file for identifyDistal3UTR
    ch_sample
        .map{ it -> it[3] }
        .combine(PREPROCESSING.out.ch_3utr_input)
        .unique()
        .set { ch_3utr_input }

    PREPROCESSING.out.ch_sample_bedgraph_files_dir
        .set{ ch_sample_bedgraph_files_dir }

    // If we are running differential, only the first inputs in the channels are required.
    // The input files used are in ch_sample_bedgraph_files_dir
    if ( run_differential ) {
        ch_3utr_input
            .first()
            .set{ ch_3utr_input }
        PREPROCESSING.out.ch_sample_bedgraph_files_dir
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

    // Convert PredictAPA output into identification and quantification challenges files
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

