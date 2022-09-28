#!/usr/bin/env python3

from __future__ import division
import io
import os
import json
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from sklearn.metrics import auc
import JSON_templates
import apaeval as apa

def main(args):

    # input parameters
    participant_input = args.participant_data
    gold_standards_dir = args.gold_standards_dir
    challenge_ids = args.challenge_ids
    participant = args.participant_name
    community = args.community_name
    out_path = args.output
    windows = args.windows
    genome_path = args.genome_dir

    print(f"INFO: input {participant_input}")
    print(f"INFO: Possible challenges {challenge_ids}")

    challenge = [c for c in challenge_ids if c.split('.')[0] in str(participant_input)][0]
    
    print(f"INFO: Selected challenge {challenge}")
    print(f"INFO: Window sizes {windows}")


    # Assuring the output path does exist
    if os.path.dirname(out_path) and not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    compute_metrics(participant_input, gold_standards_dir, challenge, participant, community, out_path, windows, genome_path)


def compute_metrics(infile, gold_standards_dir, challenge, participant, community, out_path, windows, genome_path):

    # ID prefix for assessment objects
    base_id = f"{community}:{challenge}:{participant}:"
    
    # define array that will hold the full set of assessments
    all_assessments = []

    # Dict to store metric names and corresponding variables + stderr (which is currently not computed and set to 0)
    metrics = {}

    # Vectors for AUC of precision recall curve
    p_vec = []
    r_vec = []

    # check for presence of participant input data
    assert os.path.exists(infile), "Participant file not found, please check input data."

    PD = pd.read_csv(infile, delimiter="\t", header=None, index_col=False,
        names = ['chrom_p', 'chromStart_p', 'chromEnd_p', 'name_p', 'score_p', 'strand_p'])
    # Number of PD (predicted) sites
    PD_cnt = PD.shape[0]

    # Cols for selecting unique PD sites
    PD_cols = ['chrom_p', 'chromStart_p', 'chromEnd_p', 'strand_p']

    # ground truth file
    gold_standard = os.path.join(gold_standards_dir, challenge + ".bed")
    assert os.path.exists(gold_standard), "Ground truth file not found, please check input data."

    GT = pd.read_csv(gold_standard, delimiter="\t", header=None, index_col=False,
        names = ['chrom_g', 'chromStart_g', 'chromEnd_g', 'name_g', 'score_g', 'strand_g'])
    # Number of GT sites
    GT_cnt = GT.shape[0]

    # Cols for selecting unique GT sites
    GT_cols = ['chrom_g', 'chromStart_g', 'chromEnd_g', 'strand_g']

    # genome annotation file
    genome_file = apa.select_genome_file(challenge, genome_path)
    print(f"INFO: In challenge {challenge}. Using genome file: {genome_file}")
    genome = apa.load_genome(genome_file)
    
    # For AUC of PR-curve we need more window sizes
    windowlist = list(range(0,max(windows),max(windows)//10))
    # Just in case we're missing the specified windows now
    windowlist.extend([w for w in windows if w not in windowlist])
    windowlist.sort()
    for window in windowlist:

        # Get matching sites
        # Beware: 1 PD site matching x GT sites will appear on x rows, and vice versa; Dropping duplicates is thus necessary in certain subsequent steps
        matched = apa.bedtools_window(infile, gold_standard, window)

        matched.sort_values(by=['chrom_p', 'chromStart_p', 'chromEnd_p', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])

        # METRICS for all window sizes
        ######################################

        # "true" GT: number of GT sites with matching PD sites
        tGT_cnt = matched[~matched.duplicated(GT_cols)].shape[0]
        # "true" PD: number of PD sites with matching GT sites
        tPD_cnt = matched[~matched.duplicated(PD_cols)].shape[0]

        # number of GT sites with matching PD site
        TP = tGT_cnt
        # number of prediction sites without ground truth sites
        FP = PD_cnt - tPD_cnt
        # number of ground truth sites without prediction sites
        FN = GT_cnt - TP

        sensitivity = apa.sensitivity(TP, FN)
        r_vec.append(sensitivity)
        precision = apa.precision(TP, FP)
        p_vec.append(precision)
        

        # METRICS for specified window sizes
        #########################################

        if window in windows:

            fd_rate = apa.fdr(TP, FP)
            f1_score = apa.f1_score(precision, sensitivity)
            jaccard_index = apa.jaccard(TP, FP, FN)

            # Duplicated PD sites, i.e. PD sites that matched multiple GT sites
            dPD_cnt = sum(matched.duplicated(PD_cols)) 
            dPD_proportion = dPD_cnt / tPD_cnt

            # Duplicated GT sites, i.e. GT sites that matched multiple PD sites
            dGT_cnt = sum(matched.duplicated(GT_cols))
            dGT_proportion = dGT_cnt / tGT_cnt

            # Save for assessment objects
            metrics[f"Sensitivity:{window}nt"] = [sensitivity, 0] # Metric 1
            metrics[f"Precision:{window}nt"] = [precision, 0] # Metric 2
            metrics[f"FDR:{window}nt"] = [fd_rate, 0]
            metrics[f"F1_score:{window}nt"] = [f1_score, 0]
            metrics[f"Jaccard_index:{window}nt"] = [jaccard_index, 0]
            metrics[f"multi_matched_PD:{window}nt"] = [dPD_proportion, 0] # Metric 4
            metrics[f"multi_matched_GT:{window}nt"] = [dGT_proportion, 0]

    # With precision and recall for all window sizes
    pr_auc = auc(r_vec, p_vec)

    # Save for assessment objects
    metrics["AUC"] = [pr_auc, 0] # Metric 3

    # METRICS:Percentage of genes with correctly identified number of PAS
    pas_per_gene_GT = [apa.find_pas_genes(gene, GT) for i, gene in genome.iterrows()]
    # count number of PAS per gene (i.e. row) in GT
    counts_GT = [len(np.where(pas)[0]) for pas in pas_per_gene_GT]
    # Number of genes with > 0 PAS
    genes_nonzero_PAS = np.array([x > 0 for x in counts_GT], dtype=bool)
    n_genes_nonzero_PAS = len(np.where(genes_nonzero_PAS)[0])

    # count number of PAS per gene (i.e. row) in PD
    PDwcn = PD[PD_cols].drop_duplicates()
    PDwcn.columns = GT_cols
    pas_per_gene_PD = [apa.find_pas_genes(gene, PDwcn) for i, gene in genome.iterrows()]
    counts_PD = [len(np.where(pas)[0]) for pas in pas_per_gene_PD]
    # Number of genes with identical number of polyA sites in PD and GT and > 0 PAS
    genes_identical_PAS = np.array([counts_PD[i] == counts_GT[i] for i in range(len(counts_GT))], dtype=bool)
    n_genes_identical_PAS = len(np.where(genes_identical_PAS & genes_nonzero_PAS)[0])
    perc_genes_w_PAS =  n_genes_identical_PAS / n_genes_nonzero_PAS * 100

    metrics["percentage_genes_w_correct_nPAS"] = [perc_genes_w_PAS, 0] # Metric 5

    # METRICS: polyA sites assigned to different genomic features
    # For the predicted PAS and based on the genome annotation
    exome = apa.load_genome(genome_file, feature="exon")
    prPD = apa.load_bed(infile)

    # Count terminal exons
    texome = apa.get_terminal_exons(exome)
    gtes = apa.merge_TE_overlaps(texome)
    n_tes_PD = apa.sum_overlaps(gtes, prPD)

    # Count exclusive exons, i.e. exons without terminal exons
    exclusive_exome = apa.get_non_terminal_exons(exome, texome)
    n_es_PD = apa.sum_overlaps(exclusive_exome, prPD)

    # Count introns
    introme = apa.get_introns(genome, exome)
    # report number of PAS assigned to introns
    n_intron_PD = apa.sum_overlaps(introme, prPD)

    intergenome = apa.get_intergenic_regions(gtes, genome)
    n_intergenome_PD = apa.sum_overlaps(intergenome, prPD)

    # assign remaining predictions to other
    n_other = len(prPD) - n_tes_PD - n_es_PD - n_intron_PD - n_intergenome_PD

    # for the challenge, create all assessment json objects and append them to all_assessments
    for key, value in metrics.items():
        object_id = base_id + key
        assessment_object = JSON_templates.write_assessment_dataset(object_id, community, challenge, participant, key, value[0], value[1])
        all_assessments.append(assessment_object)

    # once all assessments have been added, print to json file
    with io.open(out_path, mode='w', encoding="utf-8") as f:
        jdata = json.dumps(all_assessments,                 
                    sort_keys=True,
                    indent=4,
                    separators=(',', ': '))
        f.write(jdata)

    # print non-OEB-compatible metrics (lists) to separate file
    rogue_metrics = {
        "community_id": community,
        "challenge_id": challenge,
        "participant_id": participant}
    rogue_metrics["precision_vector"] = p_vec
    rogue_metrics["recall_vector"] = r_vec
    rogue_metrics["window_list"] = windowlist
    with io.open(os.path.join(os.path.dirname(out_path), "rogue_metrics.json"), mode='w', encoding="utf-8") as f:
        jdata = json.dumps(rogue_metrics,                 
                    sort_keys=True,
                    indent=4,
                    separators=(',', ': '))
        f.write(jdata)


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-i", "--participant_data",  help="List of execution workflow output files", required=True)
    parser.add_argument("-c", "--challenge_ids", nargs='+', help="List of challenge ids selected by the user, separated by spaces", required=True)
    parser.add_argument("-g", "--gold_standards_dir", help="dir that contains gold standard datasets for current challenge", required=True)
    parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
    parser.add_argument("-com", "--community_name", help="name/id of benchmarking community", required=True)
    parser.add_argument("-o", "--output", help="output path where assessment JSON files will be written", required=True)
    parser.add_argument("-w", "--windows", nargs='+', help="windows (nt) for scanning for poly(A) sites; several window sizes separated by spaces. The first entry is used for MSE calculation.", required=True, type=int)
    parser.add_argument("-gtf", "--genome_dir", help="genome annotation directory. Used for relative PAS usage calculation. Directory needs to contain genome files with matching organism name from challenge.", required=True)

    args = parser.parse_args()


    main(args)
