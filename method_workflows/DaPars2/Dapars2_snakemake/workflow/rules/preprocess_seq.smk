rule generateBedgraph:
    """This step generates bedgraph files from sample bams. The bedgraph files will then processed by DaPars2.
    """
    input:
        bam=lambda wildcards:
            pd.Series(samples.loc[wildcards.sample, "bam"]).values
    output:
        out=os.path.join(config["out_dir"], "{sample}.bedgraph")

    container:
        "docker://apaeval/dapars2:1.0"
    log:
        os.path.join(LOG_DIR, "{sample}_generateBedgraph.log")
    shell:
        """(genomeCoverageBed -bg -ibam {input.bam} >{output.out}) &> {log}"""


rule generatePathReadcounts:
    """This step counts the aligned read of each sample and write the abspaths and the numbers to a txt. The output will be assigned to DaPars2_config_file.txt and used for sequencing depth normalization.
    """
    input:
        bam = lambda wildcards:
            pd.Series(samples.loc[wildcards.sample, "bam"]).values
    output:
        out = os.path.join(config["out_dir"], "{sample}_pathReadcounts.tsv")
    params:
        basename = lambda wildcards: os.path.basename(samples.loc[wildcards.sample, "bam"])

    container:
        "docker://apaeval/dapars2:1.0"

    log:
        os.path.join(LOG_DIR,"{sample}_generatePathReadcounts.log")

    shell:
        """
        counts=$(samtools view -c -F 260 {input.bam})
        (echo "{params.basename}"'\t'${{counts}} > {output.out}) &> {log}
        """
