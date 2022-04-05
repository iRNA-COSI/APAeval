import os

rule download_ref_PAS:
    output:
        ref_PAS = os.path.join(config["ref_PAS_file"])

    params:
        ref_PAS_URL = config["ref_PAS_URL"]
    
    shell:
        """
        wget {params.ref_PAS_URL} -N -O {output.ref_PAS}
        """
