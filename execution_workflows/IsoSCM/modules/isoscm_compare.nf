// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/* 
    ISOSCM compare is the second step of the tool. The splice-graphs outputs from the assemble step is used in the compare
    step to perform joint segmentation of the terminal exons of coexpressed transcripts. This step reports the realtive usage
    of identified change points in each sample in a tabular format
*/
process ISOSCM_COMPARE {
        publishDir "${params.outdir}/isoscm/isoscm_compare_out_dir", mode: params.publish_dir_mode
	container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        tuple val(sample), path(assemble_out_xml_rename)

        output:
        tuple val(sample), path(compare_out), emit: ch_isoscm_compare_out

        script:
      	compare_out = sample + ".txt"
        
        """
        java -Xmx100G -jar /IsoSCM.jar compare -base $sample -x1 $assemble_out_xml_rename -x2 $assemble_out_xml_rename
        
        mv compare/$compare_out $compare_out
        """
}
