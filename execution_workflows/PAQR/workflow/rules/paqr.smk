import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import min_version
min_version("6.0")


# Download scripts required in the PAQR module
HTTP = HTTPRemoteProvider()


rule get_PAQR_scripts:
    input:
        HTTP.remote("raw.githubusercontent.com/zavolanlab/PAQR2/bd9781f1420ba39648120a6898038f49728541f3/scripts/{script}", keep_local=True)
    output:
        os.path.join(config["paqr"]["PAQ_scripts_dir"], "{script}")
    shell:
        "mv {input} {output}"


# Import the PAQR module
module paqr:
    snakefile: github("zavolanlab/PAQR2", path="Snakefile", commit="fdf56bc")
    config: config["paqr"]

use rule * from paqr as PAQR_*

use rule PAQ_create_coverages from paqr as PAQR_PAQ_create_coverages with:
    input:
        TEMP_ = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "PAQ_outdir"
        ),
        BAM_alignment = lambda wildcards:
                samples.loc[wildcards.sample_ID,"bam"],
        BAI_alignment_index = lambda wildcards:
                ".".join([samples.loc[wildcards.sample_ID,"bam"],"bai"]),
        BED_pas = os.path.join(config["out_dir"],
                str(config["atlas_version"]) + ".tpas." + config["tpas"]["strandedness"][0] + ".bed"),
        SCRIPT_ = os.path.join(
            config["paqr"]["PAQ_scripts_dir"],
            "create-pas-coverages.py"
        )

use rule PAQ_infer_relative_usage from paqr as PAQR_PAQ_infer_relative_usage with:
    input:
        BED_pas = os.path.join(config["out_dir"],
                str(config["atlas_version"]) + ".tpas." + config["tpas"]["strandedness"][0] + ".bed"),
        PKL_pas_coverage = expand(
            os.path.join(
                config["paqr"]["PAQ_outdir"],
                "pas_coverages",
                "{sample_ID}.pkl"
            ),
            sample_ID = samples.index
        ),
        TSV_extensions = expand(
            os.path.join(
                config["paqr"]["PAQ_outdir"],
                "pas_coverages",
                "{sample_ID}.extensions.tsv"
            ),
            sample_ID = samples.index
        ),
        SCRIPT_ = os.path.join(
            config["paqr"]["PAQ_scripts_dir"],
            "infer-pas-expression.py"
        )

use rule PAQ_filter_on_expression from paqr as PAQR_PAQ_filter_on_expression with:
    output:
        TSV_filtered_expression = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "filtered_pas_expression_tandem_pas.tsv"
        ),
        TSV_filtered_pas_positions = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "filtered_pas_positions_tandem_pas.tsv"
        ),
        TSV_filtered_usage = os.path.join(
            config["paqr"]["PAQ_outdir"],
            "filtered_pas_usage_tandem_pas.tsv"
        )
