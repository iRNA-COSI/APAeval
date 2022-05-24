#!/usr/bin/env nextflow

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
	   nextflow run main.nf -profile docker --input {execution.wf.APA.prediction.file} --participant_id {tool.name} --goldstandard_dir {gold.standards.dir} --challenges_ids {analyzed.challenges} --assess_dir {benchmark.data.dir} --results_dir {output.dir}
	    Mandatory arguments:
	        --input                 BED/TXT file with APA site information
	        --community_id          Name or OEB permanent ID for the benchmarking community
	        --participant_id        Name of the tool used for prediction
	        --goldstandard_dir      Dir that contains gold standard/ground truth files used to calculate the metrics for all challenges
	        --challenges_ids        Challenge ids selected by the user, separated by spaces
	        --assess_dir            Dir where the data for the benchmark are stored
	    Other options:
	        --validation_result     .json output from validation step
	        --assessment_results    .json output from metrics computation step
		--consolidation_result .json output from consolidation step
	        --outdir                The output directory where results for VRE will be saved
	        --statsdir              The output directory with nextflow statistics
	        --otherdir              The output directory where custom results will be saved (no directory inside)
	        --windows               Window sizes for scanning for poly(A) sites (List of int).
		--genome_dir            Dir where genome files for computing relative PAS usage metrics.
	        --offline               If set to 1, consolidation will be performed with local data in assess_dir only (omit to perform OEB DB query)
	    Flags:
	        --help                  Display this message
	    """.stripIndent()

	exit 1
} else {

	log.info """\

	    ==============================================
	    APAeval QUANTIFICATION BENCHMARKING PIPELINE
	    ==============================================
	        Input file: ${params.input}
	        Benchmarking community = ${params.community_id}
	        Tool name : ${params.participant_id}
	        Gold standard dataset directory: ${params.goldstandard_dir}
	        Challenge ids: ${params.challenges_ids}
	        Published benchmark data directory: ${params.assess_dir}
	        Validation result JSON file: ${params.validation_result}
	        Assessment result JSON file: ${params.assessment_results}
	        Consolidation result JSON file: ${params.consolidation_result}
	        Consolidated benchmark results directory: ${params.outdir}
	        Nextflow statistics directory: ${params.statsdir}
	        Directory with community-specific results: ${params.otherdir}
	        Window size for scanning for poly(A) sites: ${params.windows}
	        Genome dir for computing relative PAS usage metrics: ${params.genome_dir}
	        Offline mode: ${params.offline}
		""".stripIndent()

}


// input files

input_file = file(params.input)
tool_name = params.participant_id.replaceAll("\\s","_")
gold_standards_dir = Channel.fromPath(params.goldstandard_dir, type: 'dir' ) 
challenge_ids = params.challenges_ids
benchmark_data = Channel.fromPath(params.assess_dir, type: 'dir' )
community_id = params.community_id
event_date = params.event_date
windows = params.windows
genome_dir = Channel.fromPath(params.genome_dir, type: 'dir' )
genome_dir.into{
	genome_dir_val
	genome_dir_comp
}
offline = params.offline

// output 
validation_file = file(params.validation_result)
assessment_file = file(params.assessment_results)
consolidation_file = file(params.consolidation_result)
aggregation_dir = file(params.outdir, type: 'dir')
other_dir = file(params.otherdir, type: 'dir')


process validation {

	// validExitStatus 0,1
	tag "Validating input file format"
	
	publishDir "${validation_file.parent}", saveAs: { filename -> validation_file.name }, mode: 'copy'

	input:
	file input_file
	val challenge_ids
	val tool_name
	val community_id
	path genome_dir_val

	output:
	val task.exitStatus into EXIT_STAT
	file 'validation.json' into validation_out
	
	"""
	python /app/validation.py -i $input_file -com $community_id -c $challenge_ids -p $tool_name -o validation.json --genome_dir $genome_dir_val
	"""

}

process compute_metrics {

	tag "Computing benchmark metrics for submitted data"
	
	publishDir "${assessment_file.parent}", saveAs: { filename -> assessment_file.name }, mode: 'copy'

	input:
	val file_validated from EXIT_STAT
	file input_file
	val challenge_ids
	path gold_standards_dir
	val tool_name
	val community_id
	val windows
	path genome_dir_comp

	output:
	file 'assessment.json' into assessment_out

	when:
	file_validated == 0

	"""
	python3 /app/compute_metrics.py -i $input_file -c $challenge_ids -g $gold_standards_dir -p $tool_name -com $community_id -o assessment.json -w $windows --genome_dir $genome_dir_comp
	"""
}

process benchmark_consolidation {

	tag "Performing benchmark assessment and building plots"
	publishDir "${aggregation_dir.parent}", pattern: "aggregation_dir", saveAs: { filename -> aggregation_dir.name }, mode: 'copy'
	publishDir "${consolidation_file.parent}", pattern: "consolidated_result.json", saveAs: { filename -> consolidation_file.name }, mode: 'copy'

	input:
	path benchmark_data
	file assessment_out
	file validation_out
	val challenge_ids
    val event_date
	val offline
	
	output:
	path 'aggregation_dir', type: 'dir'
	path 'consolidated_result.json'

	"""
	python /app/aggregation.py -b $benchmark_data -a $assessment_out -o aggregation_dir -d $event_date --offline $offline
	python /app/merge_data_model_files.py -v $validation_out -m $assessment_out -c $challenge_ids -a aggregation_dir -o consolidated_result.json
	"""

}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
