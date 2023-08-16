// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GTFTOGENEPRED {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode
    container "docker.io/apaeval/gtftobed:1.0"

    input:
    path gtf
    
    output:
    path "*.genePred", emit: genepred
    
    script:
    genepred= 'reference.genePred'
    """
    /gtfToGenePred -genePredExt $gtf $genepred
    """
}
