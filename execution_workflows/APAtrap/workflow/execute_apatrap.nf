////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
def run_identification = modules['final_output'].run_identification
def run_quantification = modules['final_output'].run_quantification
def run_differential = modules['final_output'].run_differential

include { INPUT_CHECK } from '../subworkflows/input_check' addParams( options: [:] )
include { PREPROCESS_FILES } from '../subworkflows/preprocess_files' addParams( options: [:] )
include { RUN_IDENTIFICATION_QUANTIFICATION } from '../subworkflows/run_identification_quantification' addParams( options: [:] )
include { RUN_DIFFERENTIAL } from '../subworkflows/run_differential' addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_APATRAP{

    INPUT_CHECK ( ch_input )
        .set { ch_sample }

    PREPROCESS_FILES ( ch_sample )

    if( run_identification || run_quantification ){
        RUN_IDENTIFICATION_QUANTIFICATION (
            PREPROCESS_FILES.out.ch_sample_bedgraph_files_dir,
            PREPROCESS_FILES.out.ch_3utr_input )
    }
    if( run_differential ) {
        RUN_DIFFERENTIAL (
            PREPROCESS_FILES.out.ch_preprocessing_input,
            PREPROCESS_FILES.out.ch_sample_bedgraph_files_dir,
            PREPROCESS_FILES.out.ch_3utr_input
        )
    }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
