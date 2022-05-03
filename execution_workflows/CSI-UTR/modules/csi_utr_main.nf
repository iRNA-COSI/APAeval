// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def modules = params.modules.clone()
def csi_params = modules['csi_utr']

process CSI_UTR_MAIN {
        publishDir "${params.outdir}/csi_utr", mode: params.publish_dir_mode
        container "docker.io/apaeval/csi-utr:latest"

        input:
        tuple val(sample), val(bam)
        path sample_info
        tuple path(CSI_bed_file), path(CSI_annotation_file)

        output:
        path outdir, emit: ch_csi_utr_outdir

        script:
        genome = csi_params.genome
        r = csi_params.r
        coverage_cut = csi_params.coverage_cut
        usage_cut = csi_params.usage_cut
        p = csi_params.p
        q = csi_params.q
        bed_data_dir = "$PWD/${params.outdir}/csi_utr/input_files/"
        bed_file="./data/locations/${genome}.CSIs.bed" 
        annot_file   = "./data/annotations/${genome}.CSIs.annot.bed"
        outdir = "CSI_OUT/."
        """
        mkdir data
        mkdir data/locations
        mv $CSI_bed_file ./data/locations/.
        
        mkdir data/annotations
        mv $CSI_annotation_file ./data/annotations/.

        mkdir CSI_OUT

        CSI-UTR \
        -genome=$genome \
        -r=$r \
        -sample_info=$sample_info \
        -bed=$bed_file \
        -annot=$annot_file \
        -out=$outdir \
        -data_dir=$bed_data_dir \
        -coverage_cut=$coverage_cut \
        -usage_cut=$usage_cut \
        -p=$p \
        -q=$q \
        """
 }



