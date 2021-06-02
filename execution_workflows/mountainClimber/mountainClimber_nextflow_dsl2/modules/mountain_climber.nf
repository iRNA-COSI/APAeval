process get_junction_reads {

    container 'apaeval/exwf_mountainclimber_method:65132f2'

    input:
    path bam
    val strand
    // TODO: set strand 
      
    output:
    path './tmp/SAMPLE_junction_reads.bed', emit: junction_reads
    // TODO: set sample name (also in `script`)

    script:
    """
    get_junction_counts.py -i $bam -s $strand -o ./tmp/SAMPLE_junction_reads.bed
    """

}

process merge_tus {

    container 'apaeval/exwf_mountainclimber_method:65132f2'

    input:
    path bed
    // TODO: merge BED inputs
    path gtf
    val strand
    // TODO:
    // 1. set strand 
    // 2. there may be an issue with stranded libraries: https://github.com/gxiaolab/mountainClimber/issues/4
      
    output:
    path './tmp/transcriptional_units_merged.gtf', emit: transcriptional_units_merged

    script:
    """
    merge_tus.py -i $bed -g $gtf -s $strand -o tmp/annotations_merged
    """

}
process mountain_climber_tu {

    container 'apaeval/exwf_mountainclimber_method:65132f2'

    input:
    path bedgraph
    path junctions
    path ref_seq_sizes
    val strand
    // TODO: set strand 
      
    output:
    path './tmp/SAMPLE_transcriptional_units.bed', emit: transcriptional_units
    // TODO: set sample name (also in `script`)

    script:
    """
    mountainClimberTU.py -b $bedgraph -j $junctions -g $ref_seq_sizes -s $strand -o ./tmp/SAMPLE_transcriptional_units.bed
    """

}