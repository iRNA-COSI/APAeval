#!/usr/bin/env nextflow

nextflow.enable.dsl=2
gtf=""

process gtfToBed{
  container "docker.io/apaeval/gtftobed:1.0"
  publishDir "gtftobed"

  input:
  path gtf
    
  output:    
  path "*.bed" , emit: bed

  script:
  """
  python3 /app/gtfTobed.py --gtf $gtf
  """
}

workflow {
  gtfToBed(
    file(params.gtf)
  )
}
