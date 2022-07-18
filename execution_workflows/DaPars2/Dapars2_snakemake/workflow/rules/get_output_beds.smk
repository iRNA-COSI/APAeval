
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
        out = os.path.join(config["out_dir"], "{sample}" + config["out_id_bed_suffix"])

    params:
        out_type = "identification"

    log:
        os.path.join(LOG_DIR, "{sample}_identification_bed.log")

    container:
        "docker://amancevice/pandas:1.3.3"

    shell:
        """(python workflow/scripts/results_to_bed.py \
        -i '{input.names}' \
        -c bed \
        -t {params.out_type} \
        -o {output.out}) &> {log}"""


rule rel_usage_bed:
    """
    Take chromosome by chromosome DaPars2 results files and
    generate a BED file of identified poly(A) sites per sample
    with their fractional relative usage stored in the 'score' field
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
        out = os.path.join(config["out_dir"], "{sample}" + config["out_rel_bed_suffix"])

    params:
        out_type = "quant_relative_usage"

    log:
        os.path.join(LOG_DIR, "{sample}_rel_usage_bed.log")

    container:
        "docker://amancevice/pandas:1.3.3"

    shell:
        """(python workflow/scripts/results_to_bed.py \
        -i '{input.names}' \
        -c bed \
        -t {params.out_type} \
        -o {output.out}) &> {log}"""
