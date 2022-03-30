// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def final_output = modules['final_output']

/*
    Generate file for differential challenge
*/
process POSTPROCESS_DIFFERENTIAL {
    container "docker.io/apaeval/csi-utr:latest"

    input:
    val indicator

    output:
    path differential_out, emit: ch_differential_out

    script:
    differential_out = final_output.differential_out
    csi_utr_differential_out = "$PWD/${params.outdir}/csi_utr/CSI_OUT/DifferentialExpression/WITHIN_UTR"
    
    """
    postprocess_differential.py \
    $csi_utr_differential_out \
    $differential_out 
    """
}
