// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def inputs = modules['predict_apa']

/*
    Run the second step of APAtrap: predictAPA to infer all potential
    APA sites and estimate their corresponding usages.
*/
process PREDICT_APA {
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "docker.io/apaeval/apatrap:latest"

        input:
        tuple path(reads_bedgraph_file), path(predict_apa_input)

        output:
        path "$predict_apa_output", emit: ch_de_apa_input

        script:
        num_of_grps = inputs.g
        num_of_replicates = inputs.n
        predict_apa_output = inputs.o
        min_deg_coverage = inputs.d
        min_avg_coverage = inputs.c
        min_dist = inputs.a
        window_size = inputs.w
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
