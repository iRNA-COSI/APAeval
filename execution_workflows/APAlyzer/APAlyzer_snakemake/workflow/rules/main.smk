rule main:
    """
    A rule that runs APAlyzer.
    """

    output:
        out_main=os.path.join(config["out_dir"],'main.RData')

    params:
        outdir=config["out_dir"],

    log:
        os.path.join(LOG_DIR,"main.log")

    container:
        config["container"]

    shell:
        """(Rscript  workflow/scripts/APAlyzer_main.R \
            --dir_path {params.outdir} \
            --out_main {output.out_main]};) &> {log}"""
