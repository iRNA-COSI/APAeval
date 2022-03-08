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
def run_identification = modules['run_mode'].run_identification
def run_differential = modules['run_mode'].run_differential

include { INPUT_CHECK } from '../subworkflows/input_check' addParams( options: [:] )
include { PREPROCESS_GENOME   } from '../modules/preprocess_genome' addParams( options: [:] )
include { RUN_ISOSCM    } from '../subworkflows/run_isoscm'   addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_ISOSCM {

         INPUT_CHECK ( ch_input )
               .set { ch_sample }

         Channel
            .fromPath("${files.gtf_genome_file}")
            .set{ ch_gtf_genome_file }

        Channel
            .fromPath("${files.fasta_genome_file}")
            .set{ ch_fasta_genome_file }

        /*
            Run STAR alignment on genome files
        */
        PREPROCESS_GENOME (
            ch_gtf_genome_file,
            ch_fasta_genome_file
        )
/*
        if( run_identification ){
            RUN_IDENTIFICATION (
                PREPROCESS_FILES.out.ch_convert_to_bedgraph_out,
                PREPROCESS_FILES.out.ch_extracted_3utr_out
            )
        }

        if( run_differential ) {
            RUN_DIFFERENTIAL (
                PREPROCESS_FILES.out.ch_convert_to_bedgraph_out,
                PREPROCESS_FILES.out.ch_extracted_3utr_out
            )
        }

*/
//          RUN_ISOSCM ( ch_sample )
    }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
