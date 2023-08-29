// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def final_output = modules['final_output']

/*
    Generate file for identification challenge
*/
process POSTPROCESS_IDENTIFICATION {
    container "docker.io/apaeval/csi-utr:latest"

    input:
    tuple val(sample), val(bam), val(indicator)

    output:
    path identification_out, emit: ch_identification_out

    script:
    identification_out = sample + "_" + final_output.identification_out_suffix
    csi_utr_identification_out = "$PWD/${params.outdir}/csi_utr/CSI_OUT/COVERAGE/"
    
    """
    postprocess_identification.py \
    $csi_utr_identification_out \
    $bam \
    $identification_out 
    """
}
