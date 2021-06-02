process align_genome {

    container 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'

    input:
    path index
    path fastq
    path gtf
    // TODO: paired-end
      
    output:
    path './tmp/SAMPLE_alignments_genome', emit: alignments_genome
    // TODO:
    // 1. set sample name (also in `script`)
    // 2. check how to access file of output directory

    script:
    // TODO:
    // 1. set CPUs
    // 2. get max read length from sample table & set sjdbOverhang
    // 3. there may be an issue with the sort order: https://github.com/gxiaolab/mountainClimber/issues/1
    """
    STAR \
      --genomeDir $index \
      --readFilesIn $fastq \
      --sjdbGTFfile $gtf \
      --runThreadN 8 \
      --sjdbOverhang 99 \
      --outFileNamePrefix ./tmp/SAMPLE_alignments_genome \
      --outSAMunmapped Within \
      --outFilterType BySJout \
      --outSAMattributes NH HI AS NM MD \
      --outFilterMultimapNmax 200 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --sjdbScore 1 \
      --genomeLoad NoSharedMemory \
      --outSAMtype BAM Unsorted \
      --outSAMheaderHD @HD VN:1.4 SO:unsorted
    """

}

process align_transcriptome {

    container 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'

    input:
    path index
    path fastq
    path gtf
    // TODO: paired-end
      
    output:
    path './tmp/SAMPLE_alignments_transcriptome', emit: alignments_transcriptome
    // TODO:
    // 1. set sample name (also in `script`)
    // 2. check how to access file of output directory

    script:
    // TODO:
    // 1. set CPUs
    // 2. get max read length from sample table & set sjdbOverhang
    // 3. there may be an issue with the sort order: https://github.com/gxiaolab/mountainClimber/issues/1
    """
    STAR \
      --genomeDir $index \
      --readFilesIn $fastq \
      --sjdbGTFfile $gtf \
      --runThreadN 8 \
      --sjdbOverhang 99 \
      --outFileNamePrefix ./tmp/SAMPLE_alignments_transcriptome \
      --outFilterType BySJout \
      --outSAMattributes NH HI AS NM MD \
      --outFilterMultimapNmax 200 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --sjdbScore 1 \
      --genomeLoad NoSharedMemory \
      --outSAMtype BAM Unsorted \
      --quantMode TranscriptomeSAM \
      --outSAMheaderHD @HD VN:1.4 SO:unsorted
    """

}

process prepare_index {

    container 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'

    input:
    path fasta
      
    output:
    path './tmp/index_genome', emit: index_genome

    script: // TODO: set CPUs
    """
    STAR \
      --runMode genomeGenerate \
      --genomeFastaFiles $fasta \
      --runThreadN 6 \
      --genomeDir ./tmp/index_genome_star
    """

}
