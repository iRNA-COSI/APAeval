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
        os.path.join(config["scriptsdir"], "{script}")
    shell:
        "mv {input} {output}"


# Import the TPAS module
module tandem_pas:
    snakefile: github("zavolanlab/tandem-pas", path="Snakefile", commit="98478da")
    config: config

use rule * from tandem_pas as TPAS_*

