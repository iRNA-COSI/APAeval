// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def inputs    = modules['identify_distal_3utr']

/*
    Run the first step of APAtrap: identifyDistal3UTR to refine
    annotated 3'UTRs and identify novel 3'UTRs or 3'UTR extensions.
*/
process IDENTIFY_DISTAL_3UTR {
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "docker.io/apaeval/apatrap:latest"

        input:
        tuple path(genome_file), path(reads_bedgraph_file1), path(reads_bedgraph_file2)

        output:
        path "$utr_output", emit: ch_predictapa_input

        script:
        utr_output = "predictapa_input.bed"
        window_size = inputs.w
        extension_size = inputs.e
        min_coverage = inputs.c
        min_percentage = inputs.p
        """
        identifyDistal3UTR -i $reads_bedgraph_file1 $reads_bedgraph_file2 \
                           -m $genome_file \
                           -o $utr_output \
                           -w $window_size \
                           -e $extension_size \
                           -c $min_coverage \
                           -p $min_percentage
        """
}
