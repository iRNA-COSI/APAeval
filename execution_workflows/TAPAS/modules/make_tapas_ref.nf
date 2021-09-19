// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MAKE_TAPAS_REF{
  publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"tapas") }

  container "docker.io/apaeval/gtftobed:1.0"

  input:
  path gtf
  val out_bed
    
  output:    
  path "*.bed" , emit: bed
  path "*.txt" , emit: tapas_ref

  script:
  """
  /gtfTobed.py --gtf $gtf --out_bed $out_bed
  make_tapas_ref.py $out_bed $gtf
  """
}

