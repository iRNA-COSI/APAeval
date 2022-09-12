// Import generic module functions

//python2.7 GETUTR/getutr.py --mode makeTFfasta --gff $gff --genomefasta $fasta --lasttwoexons --librarytype RNAseq

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GETUTR_PROCESS {
        println "${params.outdir}"
        publishDir "${params.outdir}", mode: params.publish_dir_mode
        container "docker.io/apaeval/getutr:latest"
        label "process_long"

        //input:
        //tuple path(gff), path(fasta)

        //output:
        //path "*.fasta", emit: ch_tffasta
        //path "*.db", emit: ch_gffutilsdb

        script:
        """
        python2.7 /software/GETUTR/getutr.py -r /home/gregor/APAeval/data/gencode.v39.annotation.gtf -i /home/gregor/APAeval/data/SRR6795718.bam
        """
}
