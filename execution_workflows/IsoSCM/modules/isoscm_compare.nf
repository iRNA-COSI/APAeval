// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ISOSCM_COMPARE {
        publishDir "${params.outdir}/isoscm/isoscm_assemble_out_dir/${sample}", mode: params.publish_dir_mode
	container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        tuple val(sample), val(strand), path(aligned_bam), path(aligned_bai)

        output:
        tuple val(sample), path(assemble_out_identification_rename), emit: ch_isoscm_assemble_out

        script:
        assemble_out_identification_rename = "identification_out"
        
      	assemble_out_identification = sample + ".cp.gtf"
        
        """
        java -Xmx100G -jar /IsoSCM.jar assemble -bam $aligned_bam -base $sample -s $strand -dir isoscm

        mv isoscm/tmp/$assemble_out_identification $assemble_out_identification_rename
        """
}
