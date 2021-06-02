process generate_bedgraphs {

    container 'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1'

    input:
    path bam
      
    output:
    path './tmp/SAMPLE.bedgraph', emit: bedgraph
    // TODO:
    // 1. set sample name (also in `script`)
    // 2. there may be an issue with the sort order: https://github.com/gxiaolab/mountainClimber/issues/2

    script:
    """
    bedtools genomecov -trackline -bg -split -ibam $bam > /tmp/SAMPLE.bedgraph
    """

}
