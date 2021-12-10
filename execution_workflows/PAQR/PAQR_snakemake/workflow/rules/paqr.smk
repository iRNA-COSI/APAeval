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
        os.path.join(config["scriptsdir"], "{script}")
    shell:
        "mv {input} {output}"


# Import the PAQR module
module paqr:
    snakefile: github("zavolanlab/PAQR2", path="Snakefile", commit="fdf56bc")
    config: config

use rule * from paqr as PAQR_*