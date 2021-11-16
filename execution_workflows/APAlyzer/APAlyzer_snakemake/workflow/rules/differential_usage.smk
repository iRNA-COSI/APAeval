# Differential APA usage
rule differential_usage:
    """
    A rule that does the differential usage estimate of the 3' UTRs.
    """
    input:
        quants = lambda wildcards: expand(
            os.path.join(
                config["out_dir"],
                "{sample}_quant_APAlyzer.csv"),
            sample=samples.index.values,
            experiment=wildcards.experiment)
    output:
        diff_usage = os.path.join(
            config["out_dir"],
            "{experiment}_differential_usage_APAlyzer.csv")
    params:
        sample_name = "{sample}"
    log:
        os.path.join(config["local_log"], "{experiment}_differential_usage_APAlyzer.log")
    container:
        "docker://apaeval/apalyzer:latest"
    shell:
       """(Rscript  APAlyzer_differential_usage.R \
            --sample_name {params.sample_name} \
            --sample_path {input.bam} \
            --reference_genome {input.genome} \
            --out_quantification {output.out_quantification}) &> {log}"""
