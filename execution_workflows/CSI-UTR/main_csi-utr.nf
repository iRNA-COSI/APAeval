#!/usr/bin/env nextflow

// write help message here
if (params.help) {

	    log.info"""
	    ==============================================
	    APAeval pilot benchmarking pipeline
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

// define parameters (add more parameters if necessary)
if (params.sample_info_sheet) { ch_input = file(params.sample_info_sheet, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }

/*
 * Check input files indicated on the samplesheet to make sure that the input files are valid input for the downstream processes
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

def get_sample_info(LinkedHashMap info) {
    return [ info.bam_prefix, info.condition, info.sample, info.genome, info.bam_dir ]
}

// create input channels for the following process(es
 ch_samplesheet_reformat
 .splitCsv(header:true, sep:',')
 .map { get_sample_info(it) }
 .map{ it -> [ it[4], it[5] ] }
 .unique()
 .set { ch_input_the_required_ones_only }



/*
* make sample_info txt files
*/

if (!params.skip_make_sample_table){

	process make_sample_table{
		input:
		path table from params.sample_info_sheet

		output:
		  val "table_path" into ch_sample_table

		script:
		"""
		Rscript $baseDir/makeSampleTable.r $table $baseDir
		"""


	}
}


/*
* run csi-utr
*/

if (!params.skip_run_csi_utr){


	process run_csi_utr{
		input:
		val bam_path from params.input_bam_file_dir
		val x from ch_sample_table
		val r from params.r
		val p from params.p
		val q from params.q
		val anno from params.annot


		output:
		path "*" into ch_process2_outputs

		script:
		"""
		cd $baseDir/CSI-UTR/CSI-UTR_v1.1.0
      ./bin/CSI-UTR.exe\
			-out=$baseDir/results \
			-sample_info="$baseDir/sample_table.txt" \
			-annot= $anno \
			-genome="Mm10" \
			-data_dir=$bam_path \
			-r=$r \
			-p=$p \
			-q=$q

		"""


	}

}


process make_outfiles{
input:
tuple coverage from $baseDir/results/COVERAGE



}





workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
