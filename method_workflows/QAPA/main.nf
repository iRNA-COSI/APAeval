#!/usr/bin/env nextflow
/*
========================================================================================
                         PILOT BENCHMARK PIPELINE RUNNING QAPA
========================================================================================
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

if (params.help) {
    //print something here
    exit 0
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
//if (params.enable_conda) {
//    Checks.check_conda_channels(log)
//}

// Check AWS batch settings
//Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
//Checks.hostname(workflow, params, log)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    /*
     * RUN QAPA WORKFLOW
     */
    include { EXECUTE_QAPA } from './workflow/execute_qapa' addParams( [:] )
    EXECUTE_QAPA ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
