// Import generic module functions

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GETUTR_PROCESS {
        publishDir "${params.outdir}/getutr", mode: params.publish_dir_mode
        container "docker.io/apaeval/getutr:latest"
        label "process_high"

        input:
        tuple val(sample), path(bam), path(gtf)

        output:
        path "*.bed", emit: ch_getutr_output

        script:
        """
        python2.7 /software/GETUTR/getutr.py -r $gtf -i $bam -o $sample
        """
}
