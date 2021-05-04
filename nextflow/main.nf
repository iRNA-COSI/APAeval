#!/usr/bin/env nextflow

if (params.help) {
	
	    log.info"""
	    ==============================================
	    APAeval pilot benmarking pipeline
	    ==============================================
	    Usage:
	    Run the pipeline with default parameters:
	    nextflow run main.nf -profile docker --input samplesheet.csv
	    Mandatory arguments:
                --input		samplesheet
	    Flags:
                --help			Display this message
	    """.stripIndent()

	exit 1
}

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }

/*
 * Check input samplesheet
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet from ch_input
    
    output:
    path "*reformat.csv" into ch_samplesheet_reformat
    
    script:
    """
    check_samplesheet.py $samplesheet samplesheet_reformat.csv
    """
}

def get_sample_info(LinkedHashMap sample) {
    return [ sample.sample, sample.fastq1, sample.fastq2, sample.bam, sample.bai, sample.gff, sample.fasta, sample.bed, sample.mart_export ]
}

ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { get_sample_info(it) }
    .into {
        ch_sample_miso
        ch_sample_qapa
    }

/*
 * Run MISO
 */
process MISO {
    echo true
    tag "$sample"
    publishDir "${params.outdir}/miso", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(fastq1), path(fastq2), path(bam), path(bai), path(gff), path(fasta), path(bed), path(mart_export) from ch_sample_miso

    output:
    path "*" into ch_miso_outputs

    script:
    index_dir        = "indexed"
    unsummarized_dir = "unsummarized_outputs"
    """
    index_gff --index $gff $index_dir
    miso --run $index_dir $bam --output-dir $unsummarized_dir --read-len 37 --paired-end 250 15
    summarize_miso --summarize-sample $unsummarized_dir summary_output    
    """
}

/*
 * Run QAPA and SALMON
 */
process QAPA {
    echo true
    tag "$sample"
    publishDir "${params.outdir}/qapa", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(fastq1), path(fastq2), path(bam), path(bai), path(gff), path(fasta), path(bed), path(mart_export) from ch_sample_qapa

    output:
    path "*" into ch_qapa_outputs

    script:
    indexed_fasta  = "output_sequences.fa"
    utr_library    = "utr_library"
    salmon_results = "salmon_results"
    salmon_quantsf = "salmon_results/quant.sf"
    """
    qapa fasta -f $fasta $bed $indexed_fasta
    salmon index -t $indexed_fasta -i $utr_library
    salmon quant -i $utr_library -l A -1 $fastq1 -2 $fastq2 -p 4 --validateMappings --out $salmon_results
    qapa quant --db $mart_export $salmon_quantsf > pau_results.txt
    """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
