// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
     IsoSCM assmeble step is the first step in IsoSCM. It constructs a splice-graph exons,
     and segment terminal exons using the constrained change-point detection procedure.
     One of the output files contain locations of change point that can be used for
     APAeval identification challenge. Another output file is an XML file that is used for
     the second step, the compare step.
*/
process ISOSCM_ASSEMBLE {
        publishDir "${params.outdir}/isoscm/isoscm_assemble_out_dir/", mode: params.publish_dir_mode
	container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        tuple val(sample), val(strand), path(aligned_bam), path(aligned_bai)

        output:
        tuple val(sample), path(assemble_out_xml), emit: ch_isoscm_assemble_out_xml

        script:
        assemble_out_xml = sample + ".assembly_parameters.xml"
          
        """
        java -Xmx100G -jar /IsoSCM.jar assemble -bam $aligned_bam -base $sample -s $strand -dir isoscm
        
        mv isoscm/$assemble_out_xml $assemble_out_xml
        """
       
}
