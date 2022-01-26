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

include { INPUT_CHECK } from '../subworkflows/input_check' addParams( options: [:] )
include { RUN_APATRAP } from '../subworkflows/run_apatrap' addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_APATRAP{

    INPUT_CHECK ( ch_input )
        .set { ch_sample }

    RUN_APATRAP (ch_sample)

}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
