#!/usr/bin/env nextflow
/*
========================================================================================
                         BENCHMARK PIPELINE RUNNING ISOSCM
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
     * RUN ISOSCM WORKFLOW
     */
    include { EXECUTE_ISOSCM } from './workflow/execute_isoscm' addParams( [:] )
    EXECUTE_ISOSCM ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
