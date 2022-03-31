/*
    Run APAtrap to obtain differential challenge output
*/

include { IDENTIFY_DISTAL_3UTR  } from '../modules/identify_distal_3utr' addParams( options: [:] )
include { PREDICT_APA           } from '../modules/predict_apa' addParams( options: [:] )
include { CONVERT_TO_BED        } from '../modules/convert_to_bed' addParams( options: [:] )
include { DE_APA                } from '../modules/de_apa' addParams( options: [:] )
include { CONVERT_TO_TSV        } from '../modules/convert_to_tsv' addParams( options: [:] )

workflow RUN_DIFFERENTIAL {
    take:
    ch_preprocessing_input
    ch_sample_bedgraph_files_dir
    ch_3utr_input

    main:

    ch_sample_bedgraph_files_dir
        .first()
        .set{ ch_sample_bedgraph_files_dir }

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

    ch_predict_apa_input
        .first()
        .set{ ch_predict_apa_input }

    /*
        Perform second step of APAtrap: infer all potential APA sites
    */
    PREDICT_APA( ch_sample_bedgraph_files_dir,
                ch_predict_apa_input )

    /*
        Perform third step of APAtrap: detect genes having significant changes
        in APA site usage between conditions
    */
    DE_APA( PREDICT_APA.out.ch_de_apa_input )

    /*
        Convert deAPA output into differential challenge tsv file
    */
    ch_preprocessing_input
        .combine(DE_APA.out.ch_de_apa_output)
        .set{ ch_convert_to_tsv }

    CONVERT_TO_TSV( ch_convert_to_tsv )
}

