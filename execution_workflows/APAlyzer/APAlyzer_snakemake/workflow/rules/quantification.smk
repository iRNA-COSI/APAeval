# Quantification
rule quantification:
    """A rule that does the 3' UTR quantification.
    """
    input:
        bam=lambda wildcards: os.path.join(samples.loc[wildcards.sample, "bam"]),
        genome=os.path.join(config["out_dir"],"reference_genome.RData")
    output:
        out_quantification=os.path.join(config["out_dir"],"{sample}_quant_APAlyzer.csv")
    params:
        sample_name="{sample}",
        read_orientation=lambda wildcards: samples.loc[wildcards.sample, "orientation"]
    log:
        os.path.join(config["local_log"],"{sample}_{experiment}_quantification.log")
    container:
        "docker://apaeval/apalyzer:latest"

    shell:
        """(Rscript  workflow/scripts/APAlyzer_quantification.R \
             --sample_name {params.sample_name} \
             --sample_path {input.bam} \
             --reference_genome {input.genome} \
             --read_orientation {params.read_orientation} \
             --out_quantification {output.out_quantification}) &> {log}"""