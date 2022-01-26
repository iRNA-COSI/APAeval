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
        outdir = config["out_dir"],
	    read_cutoff = config["read_cutoff"],
        strandtype = config["strandtype"]

    log:
        os.path.join(LOG_DIR,"main.log")

    container:
        config["container"]

    shell:
        """(Rscript  workflow/scripts/APAlyzer_main.R \
            --dir_path {params.outdir} \
	        --read_cutoff {params.read_cutoff} \
	        --strandtype{params.strandtype} \
            --in_main {input.in_main} \
            --out_main {output.out_main};) &> {log}"""
