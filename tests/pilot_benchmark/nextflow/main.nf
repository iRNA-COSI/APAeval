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
        ch_sample_index_gff
        ch_sample_miso
        ch_sample_utr_lib
        ch_sample_qapa
    }

if (!params.skip_miso){
    ch_sample_index_gff
       .map{ it -> it[5] }
       .unique()
       .set { ch_unindexed_gff }
    /*
     * Run index_gff
     */
    process MISO_INDEX_GFF {
        publishDir "${params.outdir}/miso", mode: params.publish_dir_mode

        input:
        path gff from ch_unindexed_gff

        output:
        path "$index_dir" into ch_indexed_gff

        script:
        index_dir        = "indexed"
        """
        index_gff --index $gff $index_dir
        """
    }

    ch_indexed_gff
       .combine ( ch_sample_miso )
       .map { it -> [ it[0], it[1], it[4], it[5] ] }
       .set { ch_input_miso }

    /*
     * Run MISO
     */
    process MISO {
        tag "$sample"
        publishDir "${params.outdir}/miso/$sample", mode: params.publish_dir_mode

        input:
        tuple path(index_dir), val(sample), path(bam), path(bai) from ch_input_miso

        output:
        path "*" into ch_miso_outputs

        script:
        unsummarized_dir = "unsummarized_outputs"
        summarized_dir   = "summarized_outputs"
        """
        miso --run $index_dir $bam --output-dir $unsummarized_dir --read-len 37 --paired-end 250 15
        summarize_miso --summarize-sample $unsummarized_dir $summarized_dir    
        """
    }
}

if (!params.skip_qapa){
    ch_sample_utr_lib
       .map { it -> [ it[6], it[7] ] }
       .unique()
       .set { ch_pre_utr_lib }
    /*
     * Run QAPA & SALMON PREPARE UTRLIB
     */
    process QAPA_PREPARE_UTRLIB {
        publishDir "${params.outdir}/qapa", mode: params.publish_dir_mode

        input:
        tuple path(fasta), path(bed) from ch_pre_utr_lib

        output:
        path "$utr_library" into ch_utr_library

        script:
        indexed_fasta  = "output_sequences.fa"
        utr_library    = "utr_library"
        """
        qapa fasta -f $fasta $bed $indexed_fasta
        salmon index -t $indexed_fasta -i $utr_library
        """
    }

    ch_utr_library
       .combine ( ch_sample_qapa )
       .map { it -> [ it[0], it[1], it[2], it[3], it[9] ] }
       .set { ch_input_qapa }

    /*
     * Run QAPA & SALMON QUANT
     */
    process QAPA {
        tag "$sample"
        publishDir "${params.outdir}/qapa/$sample", mode: params.publish_dir_mode

        input:
        tuple path(utr_library), val(sample), path(fastq1), val(fastq2), path(mart_export) from ch_input_qapa

        output:
        path "*" into ch_qapa_outputs

        script:
        salmon_results = "salmon_results"
        salmon_quantsf = "salmon_results/quant.sf"
        qapa_results   = "qapa_results.txt"
        fastq_param    = ("$fastq2" == "") ? "-r $fastq1" : "-1 $fastq1 -2 $fastq2"
        """
        salmon quant -i $utr_library -l A $fastq_param -p 4 --validateMappings --out $salmon_results
        qapa quant --db $mart_export $salmon_quantsf > $qapa_results
        """
    }
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
