# Differential APA usage
rule postprocessing:
    """
    Blabla
    """

    input:
        in_postprocessing = rules.main.output.out_main

    output:
        out_postprocessing = os.path.join(config["out_dir"],'differential_challenege_output.tsv')

    log:
        os.path.join(LOG_DIR,"postprocessing.log")

    container:
        config["container"]

    shell:
       """(Rscript  workflow/scripts/APAlyzer_postprocessing.R \
            --in_postprocessing {input.in_postprocessing} \
            --out_postprocessing {output.out_postprocessing}) &> {log}"""
