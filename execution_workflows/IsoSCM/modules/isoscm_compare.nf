// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ISOSCM_COMPARE {
        publishDir "${params.outdir}/isoscm/isoscm_assemble_out_dir/${sample}", mode: params.publish_dir_mode
	container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        tuple val(sample), path(assemble_out_xml_rename)

        output:
        tuple val(sample), path(relative_usage_quantification_out), emit: ch_isoscm_compare_out

        script:
        compare_out_relative_usage_quantification_rename = "relative_usage_quantification_out"
        
      	compare_out = sample + ".txt"
        
        """
        java -Xmx100G -jar /IsoSCM.jar compare -base $sample -x1 $assemble_out_xml_rename -x2 $assemble_out_xml_rename
        
        mv compare/$compare_out $compare_out_relative_usage_quantification_rename
        """
}
