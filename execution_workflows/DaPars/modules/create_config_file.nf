// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
// get the configs for this process
def config_file_options    = modules['create_config_file']

process CREATE_CONFIG_FILE {
        publishDir "${params.outdir}/dapars", mode: params.publish_dir_mode

        output:
        path "*", emit: ch_qapa_outputs

        script:
        salmon_results = "salmon_results"
        salmon_quantsf = "salmon_results/quant.sf"
        qapa_results   = "qapa_results.txt"
        fastq_param    = ("$fastq2" == "") ? "-r $fastq1" : "-1 $fastq1 -2 $fastq2"
        """
        salmon quant -i $utr_library -l A $fastq_param -p 4 --validateMappings --out $salmon_results
        qapa quant --db $mart_export $salmon_quantsf > $qapa_results
        """
 }
