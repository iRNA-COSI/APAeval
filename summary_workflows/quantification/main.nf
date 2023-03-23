#!/usr/bin/env nextflow
nextflow.enable.dsl=2


if (params.help) {
	
	log.info"""
	    =====================================================
	    APAeval QUANTIFICATION BENCHMARKING PIPELINE
	    Author(s): Yuk Kei Wan (*1,*2), Asier Gonzalez (*3), CJ Herrmann (*4)
	    *1 iRNA COSI
	    *2 Genomic Institute of Singapore, A*STAR, Singapore
	    *3 Barcelona Supercomputing Center, Barcelone, Spain
	    *4 Biozentrum, University of Basel, Switzerland
	    =====================================================
	    Usage:
	    Run the pipeline with default parameters read from nextflow.config:
	    nextflow run main.nf -profile docker
	    Run with user parameters:
	    nextflow run main.nf -profile docker --input {execution.wf.APA.prediction.file} --participant_id {tool.name} --goldstandard_dir {gold.standards.dir} --challenges_ids {analyzed.challenges} --aggregation_dir {benchmark.data.dir} --output {output.dir}
	    Mandatory arguments:
	        --input                 List of BED/TXT files with APA site information
	        --community_id          Name or OEB permanent ID for the benchmarking community
	        --participant_id        Name of the tool used for prediction
	        --goldstandard_dir      Dir that contains gold standard/ground truth files used to calculate the metrics for all challenges
	        --challenges_ids        Challenge ids selected by the user, separated by spaces
	        --aggregation_dir            Dir where the data for the benchmark are stored
	    Other options:
	        --consolidation_result .json output from consolidation step
	        --output                The output directory where results will be saved
	        --statsdir              The output directory with nextflow statistics
	        --otherdir              The output directory where custom results will be saved (no directory inside)
	        --windows               Window sizes for scanning for poly(A) sites (List of int).
	        --genome_dir            Dir with genome files for computing relative PAS usage metrics.
	        --offline               If set to 1, consolidation will be performed with local data in aggregation_dir only (omit to perform OEB DB query)
	    Flags:
	        --help                  Display this message
	    """.stripIndent()

	exit 1
} else {

	log.info """\

	    ==============================================
	    APAeval QUANTIFICATION BENCHMARKING PIPELINE
	    ==============================================
	        Input files: ${params.input}
	        Benchmarking community = ${params.community_id}
	        Tool name : ${params.participant_id}
	        Gold standard dataset directory: ${params.goldstandard_dir}
	        Challenge ids: ${params.challenges_ids}
	        Published benchmark data directory: ${params.aggregation_dir}
	        Consolidation result JSON file: ${params.consolidation_result}
	        Consolidated benchmark results directory: ${params.output}
	        Nextflow statistics directory: ${params.statsdir}
	        Directory with community-specific results: ${params.otherdir}
	        Window size for scanning for poly(A) sites: ${params.windows}
	        Genome dir for computing relative PAS usage metrics: ${params.genome_dir}
	        Offline mode: ${params.offline}
	        TPM threshold: ${params.tpm_threshold}
		""".stripIndent()

}

// Input

input_files = Channel.fromList(params.input)
tool_name = params.participant_id.replaceAll("\\s","_")
gold_standards_dir = Channel.fromPath(params.goldstandard_dir, type: 'dir' ) 
challenge_ids = params.challenges_ids
benchmark_data = Channel.fromPath(params.aggregation_dir, type: 'dir' )
community_id = params.community_id
event_date = params.event_date
windows = params.windows
genome_dir = Channel.fromPath(params.genome_dir, type: 'dir' )
tpm_threshold = params.tpm_threshold
offline = params.offline

// Output

consolidation_file = file(params.consolidation_result)
out_dir = file(params.output, type: 'dir')
results = file(params.results, type: 'dir')
other_dir = file(params.otherdir, type: 'dir')


// Process definitions

process validation {

	// validExitStatus 0,1
	tag "Validating input file format"
	
	publishDir out_dir,
	mode: 'copy',
	overwrite: false,
	saveAs: { filename -> "validated_${input_file.baseName}.json" }


	input:
	each(path input_file)
	val challenge_ids
	val tool_name
	val community_id
	path genome_dir

	output:
	val task.exitStatus, emit: validation_status
	path validation_file, emit: validation_file
	
	"""
	python /app/validation.py -i $input_file -com $community_id -c $challenge_ids -p $tool_name -o validation_file --genome_dir $genome_dir
	"""

}

process compute_metrics {

	tag "Computing benchmark metrics for submitted data"
	
	publishDir out_dir,
	mode: 'copy',
	overwrite: false,
	pattern: "${input_file.baseName}.json",
	saveAs: { filename -> "assessments_${input_file.baseName}.json" }

	input:
	val validation_status
	each (path input_file)
	val challenge_ids
	path gold_standards_dir
	val tool_name
	val community_id
	val windows
	path genome_dir
	val tpm_threshold

	output:
	path "${input_file.baseName}.json", emit: ass_json

	when:
	validation_status == 0

	"""
	python3 /app/compute_metrics.py -i $input_file -c $challenge_ids -g $gold_standards_dir -p $tool_name -com $community_id -o "${input_file.baseName}.json" -w $windows --genome_dir $genome_dir --tpm_threshold $tpm_threshold
	"""
}

process benchmark_consolidation {

	tag "Performing benchmark assessment and building plots"

	publishDir "${results.parent}", 
	pattern: "results_dir", 
	mode: 'copy',
	overwrite: false,
	saveAs: { filename -> results.name } 

	publishDir out_dir,
	pattern: "consolidated_result.json",
	mode: 'copy',
	overwrite: false,
	saveAs: { filename -> consolidation_file.name }

	input:
	path benchmark_data
	val ass_json
	val validation_file
	val challenge_ids
    val event_date
	val offline
	
	output:
	path "results_dir"
	path "consolidated_result.json"

	"""
	python /app/aggregation.py -b $benchmark_data -a $ass_json -o results_dir -d $event_date --offline $offline
	python /app/merge_data_model_files.py -v $validation_file -m $ass_json -c $challenge_ids -a results_dir -o consolidated_result.json
	"""

}

// Workflow

workflow {
	validation(
		input_files, 
		challenge_ids, 
		tool_name, 
		community_id, 
		genome_dir
		)
	validations = validation.out.validation_file.collect()

	compute_metrics(
		validation.out.validation_status,
		input_files,
		challenge_ids,
		gold_standards_dir,
		tool_name,
		community_id,
		windows,
		genome_dir,
		tpm_threshold
		)
	assessments = compute_metrics.out.ass_json.collect()


	benchmark_consolidation(
		benchmark_data,
		assessments,
		validations,
		challenge_ids,
		event_date,
		offline
		)
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
