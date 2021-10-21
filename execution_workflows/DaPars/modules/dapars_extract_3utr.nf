// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
// get the configs for this process
def dapars_extract_3utr_options    = modules['dapars_extract_3utr']

process DAPARS_EXTRACT_3UTR {
        publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode
        container "docker.io/apaeval/dapars:latest"

        input:
        tuple path(gene_symbol_txt), path(gene_model_bed)

        output:
        path "$extracted_3utr_bed", emit: ch_extracted_3utr_output

        script:
        extracted_3utr_bed  = "final_extracted_3utr.bed"
        """
        DaPars_Extract_Anno.py -b $gene_model_bed -s $gene_symbol_txt -o $extracted_3utr_bed
        """
}
