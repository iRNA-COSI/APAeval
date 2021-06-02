process prepare_reference {

    container 'quay.io/biocontainers/rsem:1.3.3--pl5262h3198e80_2'
    // TODO: make STAR available in container; see: https://github.com/nuno-agostinho/RSEM/blob/master/Dockerfile

    input:
    path gtf
    path genome
      
    output:
    path './tmp/reference', emit: reference
    // TODO: check how to access file of output directory

    script:
    // TODO: set CPUs
    """
    rsem-prepare-reference \
      --gtf $gtf \
      -p 8 \
      --star \
      $genome \
      ./tmp/reference
    """

}

process calculate_expression {

    container 'quay.io/biocontainers/rsem:1.3.3--pl5262h3198e80_2'

    input:
    path bam
      
    output:
    path './tmp/SAMPLE_expression', emit: expression
    // TODO: set sample name (also in `script`)

    script:
    // TODO: set CPUs
    """
    rsem-calculate-expression \
      --alignments $bam \
      -p 8 \
      --paired-end \
      --append-names \
      --seed 0 \
      --estimate-rspd \
      --sampling-for-bam \
      --output-genome-bam \
      ./tmp/SAMPLE_expression
    """

}