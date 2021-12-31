/*
    Run APAtrap to obtain identification and/or quantification challenge outputs
*/

include { IDENTIFY_DISTAL_3UTR  } from '../modules/identify_distal_3utr' addParams( options: [:] )
include { PREDICT_APA           } from '../modules/predict_apa' addParams( options: [:] )
include { CONVERT_TO_BED        } from '../modules/convert_to_bed' addParams( options: [:] )

workflow RUN_IDENTIFICATION_QUANTIFICATION {
    take:
    ch_sample_bedgraph_files_dir
    ch_3utr_input

    main:
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

    /*
        Perform second step of APAtrap: infer all potential APA sites
    */
    PREDICT_APA( ch_sample_bedgraph_files_dir,
                 ch_predict_apa_input )


    // Convert PredictAPA output into identification challenges files
    CONVERT_TO_BED(PREDICT_APA.out.ch_de_apa_input)
}

