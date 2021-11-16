# Preprocessing: obtain suitable input formats
import os

rule genome_preprocessing:
    """
    A rule that creates the APA sites from a gtf file.
    """
    input:
        gtf = config["gtf"]
    output:
        out_genome = os.path.join(config["out_dir"], 'reference_genome.RData')
    # params:
    #     outdir = config["out_dir"]
    # log:
    #     os.path.join(config["local_log"], "genome_preprocessing.log")

    log:
        os.path.join(LOG_DIR, "genome_preprocessing.log")

    container:
        "docker://apaeval/apalyzer:latest"

    shell:
        """(Rscript  workflow/scripts/APAlyzer_build_reference_genome.R \
            --dir_path {params.outdir} \
            --input_gtf {input.gtf} \
            --out_reference {output.out_genome};) &> {log}"""
