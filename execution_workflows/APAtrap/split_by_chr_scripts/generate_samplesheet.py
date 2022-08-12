chromosomes = [str(x) for x in range(1,23)] + ['X', 'Y', 'M']

for chromosome in chromosomes:
    with open("P19_samplesheets/samplesheet_chr" + chromosome + ".csv", "w") as f:
        f.write(
        "condition,sample,bam,bai\n" + \
        "P19_siControl_chr" + chromosome + ",P19_siControl_R1_chr" + chromosome + ",/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918577_chr" + chromosome + ".bam,/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918577_chr" + chromosome + ".bam.bai\n" + \
        "P19_siControl_chr" + chromosome + ",P19_siControl_R2_chr" + chromosome + ",/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918578_chr" + chromosome + ".bam,/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918578_chr" + chromosome + ".bam.bai\n" + \
        "P19_siSrsf3_chr" + chromosome + ",P19_siSrsf3_R1_chr" + chromosome + ",/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918579_chr" + chromosome + ".bam,/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918579_chr" + chromosome + ".bam.bai\n" + \
        "P19_siSrsf3_chr" + chromosome + ",P19_siSrsf3_R2_chr" + chromosome + ",/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918580_chr" + chromosome + ".bam,/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918580_chr" + chromosome + ".bam.bai\n" + \
        "P19_siSrsf7_chr" + chromosome + ",P19_siSrsf7_R1_chr" + chromosome + ",/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918581_chr" + chromosome + ".bam,/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918581_chr" + chromosome + ".bam.bai\n" + \
        "P19_siSrsf7_chr" + chromosome + ",P19_siSrsf7_R2_chr" + chromosome + ",/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918582_chr" + chromosome + ".bam,/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/SRR11918582_chr" + chromosome + ".bam.bai")
