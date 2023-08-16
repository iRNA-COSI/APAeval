// Import generic module functions

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GETUTR_PROCESS {
        container "docker.io/apaeval/getutr:latest"
        label "process_high"

        input:
        tuple val(sample), path(bam), path(gtf)

        output:
        tuple val(sample), path(getutr_output), emit: ch_getutr_output

        script:
        // GETUTR outputs two file .PAVA.cps.2.0.0.bed and PAVA.smoothed.2.0.0.bed
        // We only need PAVA.cps.2.0.0.bed
        getutr_output = "${sample}.PAVA.cps.2.0.0.bed"
        """
        python2.7 /software/GETUTR/getutr.py -r $gtf -i $bam -o $sample
        """
}
