// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QAPA {
        tag "$sample"
        publishDir "${params.outdir}/qapa/$sample", mode: params.publish_dir_mode

        input:
        tuple path(utr_library), val(sample), path(fastq1), val(fastq2), path(mart_export)

        output:
        path "*", emit: ch_qapa_outputs

        script:
        salmon_results = "salmon_results"
        salmon_quantsf = "salmon_results/quant.sf"
        qapa_results   = "qapa_results.txt"
        fastq_param    = ("$fastq2" == "") ? "-r $fastq1" : "-1 $fastq1 -2 $fastq2"
        """
        /home/wanyk/Downloads/apa/tools/salmon/bin/salmon quant -i $utr_library -l A $fastq_param -p 4 --validateMappings --out $salmon_results
        qapa quant --db $mart_export $salmon_quantsf > $qapa_results
        """
 }
