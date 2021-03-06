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
def files   = modules['files']
def run_star_alignment = modules['run_mode'].run_star_alignment

include { INPUT_CHECK             } from '../subworkflows/input_check' addParams( options: [:] )
include { PREPROCESS_FILES        } from '../subworkflows/preprocess_files' addParams( options: [:] )
include { GENERATE_BAM_INDEX_FILE } from '../subworkflows/generate_bam_index_file' addParams( options: [:] )
include { RUN_ISOSCM              } from '../subworkflows/run_isoscm'   addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_ISOSCM {

        INPUT_CHECK ( ch_input )
               .set { ch_sample }
        
        if( run_star_alignment ) { 
            PREPROCESS_FILES( ch_sample )
            ch_run_isoscm_input = PREPROCESS_FILES.out.ch_aligned_bam_files
        }

        else {
            GENERATE_BAM_INDEX_FILE( ch_sample )
            ch_run_isoscm_input = GENERATE_BAM_INDEX_FILE.out.ch_bam_files
	}

        RUN_ISOSCM (
            ch_run_isoscm_input
        )       
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
