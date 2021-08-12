// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def inputs = modules['predict_apa']
def workflow_option = params.workflow.clone()
def run_differential = workflow_option['run_differential']

/*
    Run the second step of APAtrap: predictAPA to infer all potential
    APA sites and estimate their corresponding usages.
*/
process PREDICT_APA {
        tag"$sample"
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "docker.io/apaeval/apatrap:latest"

        input:
        val sample_bedgraph_files_dir
        tuple val(sample), path(reads_bedgraph_file), path(predict_apa_input)

        output:
        tuple val(sample), path(predict_apa_output), emit: ch_de_apa_input

        script:
        num_of_grps = inputs.g
        num_of_replicates = inputs.n

        min_deg_coverage = inputs.d
        min_avg_coverage = inputs.c
        min_dist = inputs.a
        window_size = inputs.w
        pwd = "$PWD/${params.outdir}/$sample_bedgraph_files_dir"
        // if run differential, all sample files have to be ran together
        if (run_differential) {
            predict_apa_output = "deapa_input.txt"
            """
            #!/bin/bash
            sample_files=""
            for file in "$pwd"/*
            do
                echo \$file
                sample_files+="\$file "
            done
            predictAPA -i \$sample_files \
                   -g $num_of_grps \
                   -n $num_of_replicates \
                   -u $predict_apa_input \
                   -o $predict_apa_output\
                   -d $min_deg_coverage \
                   -c $min_avg_coverage \
                   -a $min_dist \
                   -w $window_size
            """
        }
        // otherwise if not running differential, run sample files individually
        else {
            predict_apa_output = "deapa_input_" + sample + ".txt"
            """
            predictAPA -i $reads_bedgraph_file \
                   -g $num_of_grps \
                   -n $num_of_replicates \
                   -u $predict_apa_input \
                   -o $predict_apa_output\
                   -d $min_deg_coverage \
                   -c $min_avg_coverage \
                   -a $min_dist \
                   -w $window_size
            """
        }
}
