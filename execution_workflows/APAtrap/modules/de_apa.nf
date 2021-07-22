// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
    Run the third step of APAtrap: deAPA to detect genes having significant changes
    in APA site usage between conditions.
*/
process DE_APA {
        publishDir "${params.outdir}/apatrap", mode: params.publish_dir_mode
        container "faricazjj/apatrap"

        input:
        path de_apa_input

        output:
        path "$de_apa_output", emit: ch_de_apa_output

        script:
        de_apa_output = "de_apa_output.txt"
        """
        #!/usr/bin/Rscript

        library(deAPA)
        deAPA('$de_apa_input', '$de_apa_output', 1, 1, 1, 1, 20)
        """
}
