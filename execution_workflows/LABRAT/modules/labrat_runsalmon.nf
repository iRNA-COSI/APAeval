// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process LABRAT_RUNSALMON {
        publishDir "${params.outdir}/labrat", mode: params.publish_dir_mode
//        container "quay.io/biocontainers/labrat:0.2.2--pyhdfd78af_0"

        input:
        path maketffasta
        val all_sample
        val all_fastq1
        val all_fastq2

        output:
        path "*", emit: ch_labrat_salmon_outputs

        script:
        fastq_param    = ("$all_fastq2" == "") ? "--reads1 $all_fastq1" : "--reads1 $all_fastq1 -reads2 $all_fastq2"
        """
        LABRAT.py --mode runSalmon --librarytype RNAseq --txfasta $maketffasta $fastq_param --samplename $all_sample --threads 8
        """
}
