"""
This file contains a rule that runs APAlyzer_main.R script
using outputs from preprocessing step
"""

rule main:
    """
    A rule that runs APAlyzer.
    """
    
    input:
        in_main = rules.preprocessing.output.out_preprocessing

    output:
        out_main = os.path.join(config["out_dir"],'main.RData')

    params:
        outdir = config["out_dir"]

    log:
        os.path.join(LOG_DIR,"main.log")

    container:
        config["container"]

    shell:
        """(Rscript  workflow/scripts/APAlyzer_main.R \
            --dir_path {params.outdir} \
            --in_main {input.in_main} \
            --out_main {output.out_main};) &> {log}"""
