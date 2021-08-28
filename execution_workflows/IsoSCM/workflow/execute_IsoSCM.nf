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
if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'samplesheet not specified!' }
if (params.bamdir) { ch_bam = Channel.fromPath(params.bamdir) } else { exit 1, 'bamdir not specified!' }

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
include { RUN_ISOSCM  } from '../subworkflows/run_isoscm'  addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_ISOSCM{

         INPUT_CHECK ( ch_samplesheet )

         RUN_ISOSCM ( INPUT_CHECK.out.ch_samplesheet, ch_bam )
    }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
