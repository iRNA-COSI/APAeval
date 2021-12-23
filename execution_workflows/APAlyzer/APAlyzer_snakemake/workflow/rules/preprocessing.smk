# Preprocessing: obtain suitable input formats
import os

rule rename_gtf:
    """
    A rule that renames gtf file to the correct format for preprocessing
    """

    input:
        gtf = config["gtf"]

    output:
        gtf_renamed = os.path.join(
		config["out_dir"], 
		config["gtf_organism"]+ "." + \
		config["gtf_genome_version"] + "." + \
		config["gtf_ensemble_version"] + ".gtf")
    shell:
        "cp {input} {output}"

    shell:
        "cp {input} {output}"

rule preprocessing:
    """
    A rule that creates the APA sites from a gtf file.
    """

    input:
        gtf = rules.rename_gtf.output.gtf_renamed

    output:
        out_preprocessing = os.path.join(config["out_dir"], 'preprocessing.RData')
 
    params:
       outdir = config["out_dir"],
       sample_file = config["sample_file"]

    log:
        os.path.join(LOG_DIR, "preprocessing.log")

    container:
        config["container"]

    shell:
        """(Rscript  workflow/scripts/APAlyzer_preprocessing.R \
            --dir_path {params.outdir} \
            --sample_file_path {params.sample_file} \
            --input_gtf {input.gtf} \
            --out_preprocessing {output.out_preprocessing]};) &> {log}"""
