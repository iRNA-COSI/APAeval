# Rules to make use of the singular PAS expression from PAQR.
# first rule copied from PAQR2 repository and adjusted to work with singular_pas_expression.tsv
# These rules need to run after PAQR execution, as they depend on its scripts.

rule PAQ_normalize_expression_singular_PAS:
    """
    TPM normalize the expression values by the number of mapped reads
    """
    input:
        TSV_pas_epxression_values = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "singular_pas_expression.tsv"
        ),
        TSV_distal_sites = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "tandem_pas_expression.tsv"
        ),
        SCRIPT_ = os.path.join(
            config["paqr"]["PAQ_scripts_dir"],
            "normalize-pas-expression.py"
        )

    output:
        TSV_normalized_expression = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "singular_pas_expression_normalized.tsv"
        )

    log:
        LOG_local_stdout = os.path.join(
            config["paqr"]["PAQ_logdir"],
            "PAQ_normalize_expression_singular_PAS.stdout.log"
        ),
        LOG_local_stderr = os.path.join(
            config["paqr"]["PAQ_logdir"],
            "PAQ_normalize_expression_singular_PAS.stderr.log"
        )

    singularity:
        "docker://zavolab/mapp_base_python:1.1.1"

    shell:
        """
        python {input.SCRIPT_} \
        --tandem-pas-expression {input.TSV_pas_epxression_values} \
        --distal-pas-expression {input.TSV_distal_sites} \
        --normalized-tandem-pas-expression {output.TSV_normalized_expression} \
        1> {log.LOG_local_stdout} 2> {log.LOG_local_stderr}
        """

rule concat_paqr:
    "Concatenate tandem pas and singular pas TSV."
    input:
        TSV_filtered_expression = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "tandem_pas_expression_normalized.tsv"
        ),
        TSV_normalized_expression = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "singular_pas_expression_normalized.tsv"
        )
    output:
        TSV_filtered_expression = os.path.join(
            config["out_dir"],
            "concat_pas_expression.tsv"
        )
    singularity:
        "docker://bash:4.4.18"
    shell:
        """
        tail -n +2 {input.TSV_normalized_expression} | cat {input.TSV_filtered_expression} - > {output.TSV_filtered_expression}
        """
