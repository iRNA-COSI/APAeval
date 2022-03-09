// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ISOSCM_ASSEMBLE {
        publishDir "${params.outdir}/isoscm/", mode: params.publish_dir_mode
	container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        val aligned_bam_files_dir

        output:
        val "isoscm/aligned_bam_files", emit: ch_isoscm_assemble_out

        script:
        pwd = "$PWD/${params.outdir}/$aligned_bam_files_dir"
       
        """
        #!/bin/bash
        for folder in "$pwd"/*
        do
            filename=\$(echo \$(basename \$folder))
            sample=\$(echo "\$filename" | cut -d'.' -f1)
            strand=\$(echo "\$filename" | cut -d'.' -f2)
            bam=\$folder/\$sample.Aligned.sortedByCoord.out.bam

            java -Xmx100G -jar /IsoSCM.jar assemble -bam \$bam -base \$sample -s \$strand -dir isoscm

            mv isoscm/\$sample.assembly_parameters.xml \$folder/.
            mv isoscm/tmp/\$sample.cp.filtered.gtf \$folder/.
        done
        """
 }
