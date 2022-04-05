import os

rule paqr_tsv_to_apaeval_bed:
    input:
        TSV_filtered_expression = os.path.join(
                config["out_dir"],
                "filtered_pas_expression.tsv")

    output:
        BED_quantification = os.path.join(
                config["out_dir"],
                "{sample}_" + config["challenge_code"] + "_" + config["method"] + "_" + config["outcode"] + ".bed")

    params:
        sample = "{sample}"
    
    run:
        with open(input.TSV_filtered_expression, "r") as infile, open(output.BED_quantification, "wt") as outfile:
            # Getting sample column from header line
            header = infile.readline()
            header = [h.strip() for h in header.split("\t")]
            sample_column = header.index(params.sample)

            for line in infile:
                F = line.split("\t")
                # Getting all values for bed columns
                chromosome = F[0].strip()
                start = F[1].strip()
                stop = F[2].strip()
                ID = F[3].strip()
                tpm = F[sample_column].strip()
                strand = F[5].strip()
                # Write bed
                outfile.write("\t".join([chromosome, start, stop, ID, tpm, strand]))
                outfile.write("\n")