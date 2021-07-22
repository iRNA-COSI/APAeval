/*
    Run APAtrap
*/

include { CONVERT_TO_BEDGRAPH   } from '../modules/convert_to_bedgraph' addParams( options: [:] )
include { PREPROCESSING         } from '../modules/preprocessing' addParams( options: [:] )
include { IDENTIFY_DISTAL_3UTR  } from '../modules/identify_distal_3utr' addParams( options: [:] )
include { PREDICT_APA           } from '../modules/predict_apa' addParams( options: [:] )
include { DE_APA                } from '../modules/de_apa' addParams( options: [:] )
include { POSTPROCESSING        } from '../modules/postprocessing' addParams( options: [:] )

workflow RUN_APATRAP {
    take:
    ch_sample

    main:
    // get the bam and bai files
    ch_sample
        .map { it -> [ it[3], it[4] ] }
        .unique()
        .set { ch_preprocess }

    /*
        Preprocess the input file: convert from bam to bedgraph
    */
    CONVERT_TO_BEDGRAPH( ch_preprocess )

    /*
        Ensure that file has leading chr
    */
    PREPROCESSING( CONVERT_TO_BEDGRAPH.out.ch_bedgraph )

    // combine the genome file with the input file for identifyDistal3UTR
    ch_sample
        .map{ it -> it[7] }
        .combine(PREPROCESSING.out.ch_3utr_input)
        .set { ch_3utr_input }
    /*
        Perform first step of APAtrap: identify distal 3'utr
    */
    IDENTIFY_DISTAL_3UTR( ch_3utr_input )

    // combine the bedgraph reads with the inputfile for predictAPA
    PREPROCESSING.out.ch_3utr_input
        .combine(IDENTIFY_DISTAL_3UTR.out.ch_predictapa_input)
        .set { ch_predict_apa_input }
    /*
        Perform second step of APAtrap: infer all potential APA sites
    */
    PREDICT_APA( ch_predict_apa_input )

    /*
        Perform third step of APAtrap: detect genes having significant changes
        in APA site usage between conditions
    */
    DE_APA( PREDICT_APA.out.ch_de_apa_input )

    // combine the output files
    ch_sample
        .map{ it -> it[0] }
        .combine(DE_APA.out.ch_de_apa_output)
        .set { ch_postprocessing_input }
    /*
        Postprocess the output files
    */
    POSTPROCESSING( ch_postprocessing_input )
}

