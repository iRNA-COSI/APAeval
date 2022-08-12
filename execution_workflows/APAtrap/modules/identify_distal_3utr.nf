// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def inputs    = modules['identify_distal_3utr']
def run_differential = modules['final_output'].run_differential

/*
    Run the first step of APAtrap: identifyDistal3UTR to refine
    annotated 3'UTRs and identify novel 3'UTRs or 3'UTR extensions.
*/
process IDENTIFY_DISTAL_3UTR {

        tag "$sample"
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "docker.io/apaeval/apatrap:latest"
        label 'process_big_mem'

        input:
        val sample_bedgraph_files_dir
        tuple path(genome_file), val(sample), path(utr_input_bedgraph)

        output:
        tuple val(sample), path(utr_output), emit: ch_predictapa_input

        script:
        window_size = inputs.w
        extension_size = inputs.e
        min_coverage = inputs.c
        min_percentage = inputs.p
        pwd = "$PWD/${params.outdir}/$sample_bedgraph_files_dir"
        // if run differential, all sample files have to be ran together
        if (run_differential) {
            utr_output = "predictapa_input.bed"
            """
            #!/bin/bash
            sample_files=""
            for folder in "$pwd"/*
            do
                for file in \$folder/*
                do
                    sample_files+="\$file "
                done
            done
            identifyDistal3UTR -i \$sample_files \
                -m $genome_file \
                -o $utr_output \
                -w $window_size \
                -e $extension_size \
                -c $min_coverage \
                -p $min_percentage
            """
        }
        // otherwise if not running differential, run sample files individually
        else {
            utr_output = "predictapa_input_" + sample + ".bed"
            """
            identifyDistal3UTR -i $utr_input_bedgraph \
                    -m $genome_file \
                    -o $utr_output \
                    -w $window_size \
                    -e $extension_size \
                    -c $min_coverage \
                    -p $min_percentage
            """
        }
}
