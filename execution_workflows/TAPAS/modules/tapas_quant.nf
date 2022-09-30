// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process TAPAS_QUANT{
  publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

  container "docker.io/apaeval/tapas:1.0"

  input:
  tuple path(tapas_ref), val(sample), path(read_coverage), val(read_length)
    
  output:    
  path "refFlat*", emit: filtered_txt
  tuple val(sample), path("tapas*"), emit: quant_txt

  script:
  tapas_quant_out="tapas_quant_"+"$read_coverage"
  """
  /APA_sites_detection -ref $tapas_ref -cov $read_coverage -o $tapas_quant_out -l $read_length
  """
}

