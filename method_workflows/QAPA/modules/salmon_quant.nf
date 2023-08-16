// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SALMON_QUANT {
        tag "$sample"
        publishDir "${params.outdir}/qapa/$sample", mode: params.publish_dir_mode
        container "quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

        input:
        tuple path(utr_library), val(sample), path(fastq1), val(fastq2), path(mart_export)

        output:
        tuple val(sample), path("$salmon_quantsf"), path(mart_export), emit: ch_salmon_quant_outputs

        script:
        salmon_results = "salmon_results"
        salmon_quantsf = "salmon_results/quant.sf"
        fastq_param    = ("$fastq2" == "") ? "-r $fastq1" : "-1 $fastq1 -2 $fastq2"
        """
        salmon quant -i $utr_library -l A $fastq_param -p 4 --validateMappings --out $salmon_results
        """
 }
