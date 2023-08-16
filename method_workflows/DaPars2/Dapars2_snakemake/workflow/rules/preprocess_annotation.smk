rule getBed12:
    """
    Generate a BED12 annotation file required by DaPars2 from provided GTF transcript annotation
    """
    input:
        config["gtf"]

    output:
        os.path.join(config["out_dir"], "annotation.bed")

    container:
        "docker://apaeval/gtftobed:1.0"

    log:
        os.path.join(LOG_DIR, "getBed12.log")

    shell:
        """
        /gtfTobed.py --gtf {input} --out_bed {output}
        """


rule extractGeneIdAndName:
    """
    This step extracts gene ID and gene symbol from gtf into a two-column TSV file
    """
    input:
        gtf=config["gtf"]
    output:
        out=os.path.join(config["out_dir"], "transcriptIdSymbol.txt")
    log:
        os.path.join(LOG_DIR, "extractGeneIdAndName.log")
    shell:
        """(python workflow/scripts/extractIdSymbol.py <{input.gtf} > {output.out}) &> {log}"""


rule getRegionAnnotation:
    """This step generate 3UTR region annotation.
    """
    input:
        bed = os.path.join(config["out_dir"], "annotation.bed"),
        symbol = os.path.join(config["out_dir"], "transcriptIdSymbol.txt")

    output:
        out = os.path.join(config["out_dir"], "extracted_3UTR.bed")

    container:
        "docker://apaeval/dapars2:1.0"

    log:
        os.path.join(LOG_DIR, "getRegionAnnotation.log")

    shell:
        """(python /DaPars2/src/DaPars_Extract_Anno.py \
        -b {input.bed} \
        -s {input.symbol} \
        -o {output.out}) &> {log}
        """
