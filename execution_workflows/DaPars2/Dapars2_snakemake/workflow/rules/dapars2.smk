rule makeDapars2Config:
    """
    This step makes configure file for DaPars2.
    """
    input:
        bed = os.path.join(config["out_dir"], "extracted_3UTR.bed"),
        seqDepth = os.path.join(config["out_dir"], "{sample}_pathReadcounts.tsv"),
        bedgraph = os.path.join(config["out_dir"], "{sample}.bedgraph")
    output:
        out=os.path.join(config["out_dir"], "{sample}_Dapars2_configure_file.txt")
    params:
        ct=config["coverageThreshold"],
        bedgraph = os.path.abspath(os.path.join(config["out_dir"], "{sample}.bedgraph")),
        seqDepth = os.path.abspath(os.path.join(config["out_dir"], "{sample}_pathReadcounts.tsv")),
        bed = os.path.abspath(os.path.join(config["out_dir"],  "extracted_3UTR.bed")),
        out = os.path.join(config["out_dir"], "intermediate_{sample}", "apa"),
        outfile = "apa"
    threads: 16
    log:
        os.path.join(LOG_DIR, "{sample}_makeDapars2Config.log")
    run:
        ct=params.ct
        thr=threads
        seqd = params.seqDepth
        with open(output.out, "w") as f:
            f.write("# The following file is the result of generate_region_annotation\n")
            f.write(f"Annotated_3UTR={str(params.bed)}\n")
            f.write(f"#A comma-separated list of bedgraph files of all samples\n")
            f.write(f"Aligned_Wig_files={params.bedgraph}\n")
            f.write(f"Output_directory={params.out}\n")
            f.write(f"Output_result_file={params.outfile}\n")
            f.write(f"Coverage_threshold={ct}\n")
            f.write(f"Num_Threads = {thr}\n")
            f.write(f"sequencing_depth_file={seqd}\n")



def getchromosomes(filename):
    chrs = set()
    with open(filename, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            else:
                chrs.add(line.split('\t')[0])
    chr_names = list(chrs)
    return chr_names

chromosomes = getchromosomes(config['gtf'])


rule mainDapars2:
    """
    Run DaPars2 on per chromosome for each sample
    """
    input:
        config = os.path.join(config["out_dir"], "{sample}_Dapars2_configure_file.txt")

    output:
        os.path.join(config["out_dir"], "intermediate_{sample}", "apa_{chr}", "apa_result_temp.{chr}.txt")

    params:
        chromosome = "{chr}"

    wildcard_constraints:
        chr = "|".join(chromosomes)

    container:
        "docker://apaeval/dapars2:1.0"

    log:
        os.path.join(LOG_DIR, "{sample}_{chr}_execute.mainDapars2.log")
    shell:
        """
        python /DaPars2/src/Dapars2_Multi_Sample.py \
        {input.config} \
        {params.chromosome} \
        &> {log}
        """
