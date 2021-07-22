// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Run the first step of APAtrap: identifyDistal3UTR to refine
    annotated 3'UTRs and identify novel 3'UTRs or 3'UTR extensions.
*/
process IDENTIFY_DISTAL_3UTR {
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "docker.io/apaeval/apatrap:latest"

        input:
        tuple path(genome_file), path(reads_bedgraph_file)

        output:
        path "$utr_output", emit: ch_predictapa_input

        script:
        utr_output = "predictapa_input.bed"
        """
        identifyDistal3UTR -i $reads_bedgraph_file -m $genome_file -o $utr_output
        """
}
