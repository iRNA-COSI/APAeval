#!/usr/bin/env nextflow

nextflow.enable.dsl=2
gtf=""
out_bed=""

process gtfToBed{
  container "docker.io/apaeval/gtftobed:1.0"
  publishDir "gtftobed"

  input:
  path gtf
  val out_bed
    
  output:    
  path "*.bed" , emit: bed

  script:
  """
  /gtfTobed.py --gtf $gtf --out_bed $out_bed
  """
}

workflow {
  gtfToBed(
    file(params.gtf), params.out_bed
  )
}
