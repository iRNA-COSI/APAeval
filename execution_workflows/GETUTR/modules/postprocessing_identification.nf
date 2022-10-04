// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
// get the configs for this process
def options = modules['final_output']

/*
    Convert GETUTR output file to differential challenge file
*/
process POSTPROCESSING_IDENTIFICATION {
    publishDir "${params.outdir}/getutr", mode: params.publish_dir_mode
    container "docker.io/apaeval/getutr:latest"

    input:
    val sample
    path "*"

    script:
    file1 = "$PWD/${params.outdir}/getutr/${sample}.PAVA.cps.2.0.0.bed"
    file2 = "$PWD/${params.outdir}/getutr/${sample}.PAVA.smoothed.2.0.0.bed" // dont need this one
    identification_output = "$PWD/${params.outdir}/getutr/${sample}_identification_output.bed" // dont need this one
    """
    awk '{if (\$1=="track") print \$0; else print \$1, \$2, \$3, \$4, ".", \$6}' ${file1} > temp1
    mv temp1 ${identification_output}
    rm ${file1}
    rm ${file2}
    """
}