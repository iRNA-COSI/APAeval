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

rule preprocessing:
    """
    A rule that creates the APA sites from a gtf file.
    """

    input:
        gtf = gtf_renamed

    output:
        out_preprocessing = os.path.join(config["out_dir"], 'preprocessing.RData')
 
    params:
       outdir = config["out_dir"],
       sample_file = config["sample_file"]

    log:
        os.path.join(LOG_DIR, "preprocessing.log")

    container:
        "docker://apaeval/apalyzer:1.0.3"

    shell:
        """(Rscript  workflow/scripts/APAlyzer_preprocessing.R \
            --dir_path {params.outdir} \
            --sample_file_path {params.sample_file} \
            --input_gtf {input.gtf} \
            --out_preprocessing {output.out_preprocessing]};) &> {log}"""
