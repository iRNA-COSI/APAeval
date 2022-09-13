// Import generic module functions

//python2.7 GETUTR/getutr.py --mode makeTFfasta --gff $gff --genomefasta $fasta --lasttwoexons --librarytype RNAseq

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GETUTR_PROCESS {
        publishDir "${params.outdir}", mode: params.publish_dir_mode
        container "docker.io/apaeval/getutr:latest"
        label "process_high"

        input:
        tuple val(sample), path(bam), path(gtf)

        output:
        path "*.bed"

        script:
        """
        python2.7 /software/GETUTR/getutr.py -r $gtf -i $bam -o $sample
        """
}
