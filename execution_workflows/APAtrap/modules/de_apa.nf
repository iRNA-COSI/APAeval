// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def inputs    = modules['de_apa']

/*
    Run the third step of APAtrap: deAPA to detect genes having significant changes
    in APA site usage between conditions.
*/
process DE_APA {
        tag "$sample"
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "docker.io/apaeval/apatrap:latest"
        label 'process_high'

        input:
        tuple val(sample), path(de_apa_input)

        output:
        tuple val(sample), path(de_apa_output), emit: ch_de_apa_output

        script:
        de_apa_output = "depa_output.txt"
        group1 = inputs.group1
        group2 = inputs.group2
        least_qualified_num_in_group1 = inputs.least_qualified_num_in_group1
        least_qualified_num_in_group2 = inputs.least_qualified_num_in_group2
        coverage_cutoff = inputs.coverage_cutoff
        """
        #!/usr/bin/Rscript

        library(deAPA)
        deAPA('$de_apa_input', '$de_apa_output', $group1, $group2, $least_qualified_num_in_group1, $least_qualified_num_in_group2, $coverage_cutoff)
        """
}
