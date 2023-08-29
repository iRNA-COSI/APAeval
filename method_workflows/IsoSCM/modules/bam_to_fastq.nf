// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Convert input BAM file to fastq for star alignment step
*/
process BAM_TO_FASTQ {
    container "docker.io/apaeval/isoscm:latest"

    input:
    tuple val(sample), path(bam), val(read_type)

    output:
    tuple val(sample), path(compressed_fastq1), path(compressed_fastq2), emit: ch_fastq_files

    script:
    sorted_bam = sample + ".qsort.bam"
    fastq1 = sample + ".fastq1.fastq"
    fastq2 = sample + ".fastq2.fastq"
    compressed_fastq1 = fastq1 + ".gz"
    compressed_fastq2 = fastq2 + ".gz"
    if (read_type == "paired") {
    	"""
    	samtools sort -n -o $sorted_bam $bam
        bedtools bamtofastq -i $sorted_bam \
                            -fq $fastq1 \
                            -fq2 $fastq2
        gzip -f $fastq1
        gzip -f $fastq2
    	"""
    }
    else {
        """
        bedtools bamtofastq -i $bam -fq $fastq1
        touch $compressed_fastq2
        gzip -f $fastq1
        """
    }
}
