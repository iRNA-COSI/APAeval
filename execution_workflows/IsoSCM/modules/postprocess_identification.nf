// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
def modules = params.modules.clone()
def output_files = modules['output_files'].clone()

process POSTPROCESS_IDENTIFICATION {
        publishDir "${params.outdir}/isoscm/", mode: params.publish_dir_mode
        container "docker.io/apaeval/isoscm:latest"
        label 'process_high'

        input:
        val aligned_bam_files_dir

        output:
        val "isoscm/aligned_bam_files", emit: ch_isoscm_assemble_out

        script:
        output_dir = "${output_files.output_dir}"
        pwd = "$PWD/${params.outdir}/$aligned_bam_files_dir"
        output_folder = "$PWD/${params.outdir}/$output_dir/"
        identification_out_suffix = "${output_files.identfication_out_suffix}"
        
        """
        #!/bin/bash
        for folder in "$pwd"/*
        do
            filename=\$(echo \$(basename \$folder))
            sample=\$(echo "\$filename" | cut -d'.' -f1)
            strand=\$(echo "\$filename" | cut -d'.' -f2)
            identification_file=\$folder/\$sample.cp.filtered.gtf
            output_file=\${output_folder}\${sample}_\${identification_out_suffix}
    
            postprocess_identification.py \$identification_file \$output_file
        done
        """
 }
