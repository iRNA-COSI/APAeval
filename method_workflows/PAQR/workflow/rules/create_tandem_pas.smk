import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import min_version
min_version("6.0")


# Download scripts required in the TPAS module
HTTP = HTTPRemoteProvider()

rule get_TPAS_scripts:
    input:
        HTTP.remote("raw.githubusercontent.com/zavolanlab/tandem-pas/bd46f60ef77faa46caf511ac407824bade4236c6/scripts/{script}", keep_local=True)
    output:
        os.path.join(config["tpas"]["scriptsdir"], "{script}")
    shell:
        "mv {input} {output}"


# Import the TPAS module
module tandem_pas:
    snakefile: github("zavolanlab/tandem-pas", path="Snakefile", commit="7f12eb2")
    config: config["tpas"]

use rule * from tandem_pas as TPAS_*

use rule select_tandem_pas from tandem_pas as TPAS_select_tandem_pas with:
    input:
        BED_pas_atlas = config["ref_PAS_file"],
        GTF_annotation = config["tpas"]["gtf"],
        SCRIPT_ = os.path.join(
            config["tpas"]['scriptsdir'],
            "mz-select-pas-subset.pl")

use rule filter_on_ambiguous_annotation from tandem_pas as TPAS_filter_on_ambiguous_annotation with:
    output:
        BED_tandem_pas_terminal_exons_clean = os.path.join(config["out_dir"],
                str(config["atlas_version"]) + ".tpas.{strandedness}.bed")
