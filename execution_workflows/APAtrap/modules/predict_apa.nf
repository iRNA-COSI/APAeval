// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Run the second step of APAtrap: predictAPA to infer all potential
    APA sites and estimate their corresponding usages.
*/
process PREDICT_APA {
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "faricazjj/apatrap"

        input:
        tuple path(reads_bedgraph_file), path(predict_apa_input)

        output:
        path "$predict_apa_output", emit: ch_de_apa_input

        script:
        predict_apa_output = "de_apa_input.txt"
        """
        predictAPA -i $reads_bedgraph_file -g 1 -n 1 -u $predict_apa_input -o $predict_apa_output
        """
}
