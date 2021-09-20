
rule identification_bed:
    """
    Take chromosome by chromosome DaPars2 results files and
    generate a BED file of identified poly(A) sites per sample
    """
    input:
        names = lambda wildcards: expand(os.path.join(
                            config["out_dir"],
                            "intermediate_{sample}",
                            "apa_{chr}",
                            "apa_result_temp.{chr}.txt"),
                        sample=wildcards.sample,
                        chr=chromosomes)
    output:
        out = os.path.join(config["out_dir"], "{sample}" + config["out_bed_suffix"])

    log:
        os.path.join(LOG_DIR, "{sample}_identification_bed.log")

    shell:
        """(python workflow/scripts/results_to_bed.py \
        -i '{input.names}' \
        -c bed \
        -o {output.out}) &> {log}"""
