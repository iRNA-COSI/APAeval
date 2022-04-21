import os

rule paqr_tsv_to_apaeval_bed:
    input:
        TSV_pas_expression = os.path.join(
                config["out_dir"],
                "concat_pas_expression.tsv")

    output:
        BED_quantification = os.path.join(
                config["out_dir"],
                "{sample}_" + config["challenge_code"] + "_" + config["method"] + "_" + config["outcode"] + ".bed")

    params:
        sample = "{sample}"
    
    run:
        with open(input.TSV_pas_expression, "r") as infile, open(output.BED_quantification, "wt") as outfile:
            # Getting sample column from header line
            header = infile.readline()
            header = [h.strip() for h in header.split("\t")]
            sample_column = header.index(params.sample)

            for line in infile:
                F = line.split("\t")

                # Getting all values for BED columns

                # From PAS ID (chr:rep:strand) for single nucleotide PAS
                ID = F[3].strip().split(":")
                chromosome = ID[0]
                # convert from ID's 1- to BED's 0-based coordinates
                start = str(int(ID[1]) - 1) 
                stop = ID[1] # = start + 1
                strand = ID[2]

                # From column for current sample
                tpm = F[sample_column].strip()
                
                # Write bed
                if float(tpm) != float(-1.0):
                    outfile.write("\t".join([chromosome, start, stop, ":".join(ID), tpm, strand]))
                    outfile.write("\n")
