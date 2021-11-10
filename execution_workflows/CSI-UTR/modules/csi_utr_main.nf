// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
def csi_params = modules['csi_utr']

process CSI_UTR_MAIN {
        publishDir "${params.outdir}/csi_utr", mode: params.publish_dir_mode
        container "docker.io/apaeval/csi-utr:latest"

        input:
        path sample_info
        path genome_file

        output:
        path config_file, emit: ch_dapars_output

        script:
        genome = csi_params.genome
        r = csi_params.r
        coverage_cut = csi_params.coverage_cut
        usage_cut = csi_params.usage_cut
        p = csi_params.p
        q = csi_params.q
        outdir = "$PWD/${params.outdir}/csi_utr/"
        bed_data_dir = "$PWD/${params.outdir}/csi_utr/input_files/"
        """
        CSI-UTR \
        -genome=$genome \
        -r=$r \
        -sample_info=$sample_info \
        -out=$outdir \
        -data_dir=$bed_data_dir \
        -coverage_cut=$coverage_cut \
        -usage_cut=$usage_cut \
        -p=$p \
        -q=$q \
        """
 }



