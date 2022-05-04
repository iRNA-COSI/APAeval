// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet
    
    output:
    path "*reformat.csv"
    
    script:
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}

def get_sample_info(LinkedHashMap sample) {
    return [ sample.sample, sample.bam, sample.strand, sample.read_type ]
}
