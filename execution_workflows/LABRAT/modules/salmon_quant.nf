// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SALMON_QUANT {
        tag "$sample"
        publishDir "${params.outdir}/labrat/salmon/", mode: params.publish_dir_mode
        container "quay.io/biocontainers/salmon:1.6.0--h84f40af_0"

        cpus 8

        input:
        tuple path(txfasta_idx), val(sample), path(fastq1), val(fastq2)

        output:
        tuple val(sample), path("$sample"), emit: ch_salmon_quant_outputs
        tuple val(sample), val("labrat/salmon/"), emit: ch_salmon_dir

        script:
        if ("$fastq2" == ""){
            """
            salmon quant --libType A -p ${task.cpus} --fldMean 250 --fldSD 20 --seqBias -r $fastq1 -o $sample --index $txfasta_idx --validateMappings
            """
        } else {
            """
            salmon quant --libType A -p ${task.cpus} --seqBias --gcBias -1 $fastq1 -2 $fastq2 -o $sample --index $txfasta_idx --validateMappings
            """
        }
}
