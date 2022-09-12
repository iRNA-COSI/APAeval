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
//if (params.labrat_quant_dir) { ch_labrat_quant_dir = file(params.labrat_quant_dir) } else { ch_labrat_quant_dir = '' }

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

// TO DO -- define options for the processes below

include { INPUT_CHECK } from '../subworkflows/input_check' addParams( options: [:] )
include { RUN_GETUTR  } from '../subworkflows/run_getutr'  addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_GETUTR{
         INPUT_CHECK ( ch_input )
         //RUN_GETUTR ( INPUT_CHECK.out.ch_sample, INPUT_CHECK.out.ch_condition, ch_labrat_quant_dir)
         RUN_GETUTR ( INPUT_CHECK.out.ch_sample )
    }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
