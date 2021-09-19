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

// TO DO -- define options for the processes below
def samtools_options = modules['samtools']
def tapas_quant_options    = modules['tapas_quant']

include { INPUT_CHECK    } from '../subworkflows/input_check' addParams( options: [:] )
include { SAMTOOLS       } from '../modules/samtools'         addParams( options: samtools_options )
include { MAKE_TAPAS_REF } from '../modules/make_tapas_ref'   addParams( options: [:] )
include { TAPAS_QUANT    } from '../modules/tapas_quant'      addParams( options: tapas_quant_options )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow EXECUTE_TAPAS{

         INPUT_CHECK ( ch_input )
               .set { ch_sample }

         ch_sample
               .map { it -> it[2] }
               .unique()
               .set { ch_gtf }
         MAKE_TAPAS_REF ( ch_gtf, "tapas.bed" )

         SAMTOOLS ( ch_sample )

         MAKE_TAPAS_REF.out.tapas_ref
            .combine ( SAMTOOLS.out.read_coverage )
            .set { ch_tapas_input }
         TAPAS_QUANT( ch_tapas_input )
    }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
