// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process

/*
    Convert GETUTR output file to identification challenge file
*/
process POSTPROCESSING_IDENTIFICATION {
    publishDir "${params.outdir}/getutr", mode: params.publish_dir_mode
    container "docker.io/apaeval/getutr:latest"

    input:
    tuple val(sample), path(getutr_output)

    output:
    path "*"

    script:
    identification_output = "${sample}_identification_output.bed"
    """
    awk 'BEGIN {OFS="\t"} {if (\$1=="track") print \$0; else print \$1,\$2,\$3,\$4,".",\$6}' ${getutr_output} > ${identification_output}
    """
}
