"""Snakemake pipeline for MISO.

Used for APAeval challenge
"""

import pandas as pd

samples = pd.read_csv(os.path.abspath(config["samples"])).set_index("sample", drop=False)

rule finish:
    """Rule that specifies the final output.
        In this case a json file with compute resources used (OUTCODE 04)
    """
    input:
        os.path.join(config["out_dir"],config["benchmarks"], config["param_code"], "_".join(config["param_code"], config["method"], "04.json")

###########
# Preprocessing: obtain suitable input formats

###########
# Method-specific rules

rule index:
    input:
        config["genome"]
    output:
        os.path.join(config["out_dir"], "indexed", "genes.gff")
    params:
        dir=os.path.join(config["out_dir"], "indexed")
    conda: 
        os.path.join(config["envs"], "miso.yaml")
    benchmark:
        os.path.join(config["out_dir"], config["benchmarks"], "index.tsv")
    log:
        os.path.join(config["out_dir"], config["log_dir"], "index.log")
    shell:
        "index_gff --index {input} {params.dir} &> {log}"

# TODO: use 'strandedness' to execute proper miso run
rule execute:
    input:
        index=os.path.join(config["out_dir"], "indexed", "genes.gff"),
        bam=lambda wildcards:
                pd.Series(
                    samples.loc[wildcards.sample, "bam"]
                ).values
    output:
        os.path.join(config["out_dir"], "{sample}", "summary", "{sample}.miso_summary")
    params:
        index_dir=os.path.join(config["out_dir"], "indexed"),
        settings=config["miso_settings"],
        out_dir=os.path.join(config["out_dir"], "{sample}"),
        read_len=config["read_len"]
    threads: 4
    conda: 
        os.path.join(config["envs"], "miso.yaml")
    benchmark:
        os.path.join(config["out_dir"], config["benchmarks"], "execute.{sample}.tsv")
    log:
        run=os.path.join(config["out_dir"], config["log_dir"], "run.{sample}.log"),
        summarize=os.path.join(config["out_dir"], config["log_dir"], "summarize.{sample}.log")
    shell:
        """miso --run {params.index_dir} {input.bam} \
            --settings-filename {params.settings} \
            --output-dir {params.out_dir} \
            --read-len {params.read_len} &> {log.run}; \
        summarize_miso --summarize-sample \
            {params.out_dir} {params.out_dir} &> {log.summarize}
        """

################
# Postprocessing & benchmark gathering

rule gather_benchmark_Q1:
    """Obtain runtime and max memory usage

    Per sample, obtain the runtime and max memory from the benchmarked file
    and compute sum of individual runtimes and max of all max memories.
    
    """
    input:
        T1=os.path.join(config["out_dir"], config["benchmarks"], "index.tsv"),
        T2=expand(os.path.join(config["out_dir"], config["benchmarks"], "execute.{sample}.tsv"),
            sample = samples.index)
    output:
        os.path.join(config["out_dir"],config["benchmarks"], config["param_code"], "_".join(config["param_code"], config["method"], "04.json")
    run:
        import pandas as pd
        import json
        res = {'run_time_sec': 0, 'max_mem_mib': 0}
        for file in input:
            df = pd.read_table(file, sep="\t", header = 0)
            res['run_time_sec'] += df.s.values.mean()
            max_mem = df.max_pss.max()
            if type(max_mem) is not int:
                continue
            if max_mem > res['max_mem_mib']:
                res['max_mem_mib'] = max_mem 
        with open(output[0], 'w') as json_file:
            json.dump(res, json_file, indent = 4)

