# Differential APA usage
rule postprocessing:
    """
    Blabla
    """

    output:
        out_final = os.path.join(config["out_dir"],'differential_challenege_output.tsv')

    log:
        os.path.join(LOG_DIR,"postprocessing.log")

    container:
        "docker://apaeval/apalyzer:1.0.3"

    shell:
       """(Rscript  APAlyzer_postprocessing.R \
            --out_final {output.out_final}) &> {log}"""
