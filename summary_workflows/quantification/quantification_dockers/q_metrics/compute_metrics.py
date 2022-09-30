#!/usr/bin/env python3

from __future__ import division
import io
import os
import json
import pandas as pd
from argparse import ArgumentParser
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
    tpm_threshold = args.tpm_threshold

    print(f"INFO: input {participant_input}")
    print(f"INFO: Possible challenges {challenge_ids}")

    challenge = [c for c in challenge_ids if c.split('.')[0] in str(participant_input)][0]
    
    print(f"INFO: Selected challenge {challenge}")


    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    compute_metrics(participant_input, gold_standards_dir, challenge, participant, community, out_path, windows, genome_path, tpm_threshold)


def compute_metrics(infile, gold_standards_dir, challenge, participant,
    community, out_path, windows, genome_path, TPM_THRESHOLD):

    # define array that will hold the full set of assessment datasets
    all_assessments = []
    # define list with keywords for return type of dataframe from match_with_gt()
    all_return_df_types = ["all_GT", "union", "intersection"]

    # ID prefix for assessment objects
    base_id = f"{community}:{challenge}:{participant}:"
    # Dict to store metric names and corresponding variables + stderr (which is currently not computed and set to 0)
    metrics = {}

    # check for presence of participant input data
    assert os.path.exists(infile), "Participant file not found, please check input data."

    # ground truth file
    gold_standard = os.path.join(gold_standards_dir, challenge + ".bed")
    assert os.path.exists(gold_standard), "Ground truth file not found, please check input data."

    # genome annotation file
    genome_file = apa.select_genome_file(challenge, genome_path)
    # log to stdout
    print(f"INFO: In challenge {challenge}. Using genome file: {genome_file}")
    genome = apa.load_genome(genome_file)

    # Cols for selecting unique PD sites
    PD_cols = ['chrom_p', 'chromStart_p', 'chromEnd_p', 'strand_p']

    # Cols for selecting unique GT sites
    GT_cols = ['chrom_g', 'chromStart_g', 'chromEnd_g', 'strand_g']

    for window in windows:

        ## Get matching sites
        matched = apa.bedtools_window(infile, gold_standard, window)

        # filter TPM <= TPM_THRESHOLD
        matched = matched.loc[matched.score_p > TPM_THRESHOLD,]

        # Number of duplicated PD sites, i.e. PD sites that matched multiple GT sites
        dPD_count = sum(matched.duplicated(PD_cols, keep="first")) 

        # Number of "duplicated" GT sites, i.e. GT sites that matched multiple PD sites
        dGT_count = sum(matched.duplicated(GT_cols, keep="first")) 
    
        # split PD sites that matched with multiple GT
        matched = apa.split_pd_by_dist(matched)
        if matched.empty:
            raise RuntimeError(f"No overlap found between participant: {infile} and ground truth: {gold_standard}")

        # sort
        matched.sort_values(by=['chrom_p', 'chromStart_p', 'chromEnd_p', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])

        # merge PD sites that matched with the same GT by summing their expression
        matched = apa.merge_pd_by_gt(matched)

        ## Get PD sites without matching GT (FP)
        only_PD = apa.bedtools_window(infile, gold_standard, window, reverse=True)
        # filter TPM <= TPM_THRESHOLD
        only_PD = only_PD.loc[only_PD.score_p > TPM_THRESHOLD,]
        if not only_PD.empty:
            # Duplicate the cols to cols with label "g"
            only_PD[['chrom_g', 'chromStart_g', 'chromEnd_g', 'name_g', 'score_g', 'strand_g']] = only_PD[['chrom_p', 'chromStart_p', 'chromEnd_p', 'name_p', 'score_p', 'strand_p']]
            # Now GT and PD cols are the same: FP; so label "g" has to get expression zero
            only_PD['score_g'] = [0.0]*len(only_PD)

        ## Get GT sites without matching PD (FN)
        # Note that columns at first will be labelled "p", although they are ground truth sites
        only_GT = apa.bedtools_window(gold_standard, infile, window, reverse=True)
        if not only_GT.empty:
            # Duplicate the cols to cols with label "g"
            only_GT[['chrom_g', 'chromStart_g', 'chromEnd_g', 'name_g', 'score_g', 'strand_g']] = only_GT[['chrom_p', 'chromStart_p', 'chromEnd_p', 'name_p', 'score_p', 'strand_p']]
            # Now GT and PD cols are the same: FN; so label "p" has to get expression zero
            only_GT['score_p'] = [0.0]*len(only_GT)


        # METRIC: Expression unmatched sites (= expression FP)
        ####################################
        # Key: exact name of metric as it appears in specification
        metric_name = f"Sum_FP_TPM:{window}nt"
        
        # Get metric
        if not only_PD.empty:
            # total expression of non-matched PD sites
            tpm_unmatched = (only_PD["score_p"]).sum()
        else:
            tpm_unmatched = 0.0
        
        # Value: List of [variable_holding_metric, std_err]
        metrics[metric_name] = [tpm_unmatched, 0]
        
        
        # METRIC: Percent expression unmatched
        ######################################
        # Key: exact name of metric as it appears in specification
        metric_name = f"Percent_FP_TPM:{window}nt"
        
        # Get metric
        if not matched.empty:
            tpm_matched = (matched["score_p"]).sum()
            pct_tpm_unmatched = tpm_unmatched / (tpm_unmatched + tpm_matched) * 100
        else:
            pct_tpm_unmatched = 100.0

        # Value: List of [variable_holding_metric, std_err]
        metrics[metric_name] = [pct_tpm_unmatched, 0]

        # METRIC Performance metrics
        ######################################
        # binary classification metrics
        # number of ground truth sites with matching prediction site
        TP = matched[~matched.duplicated(GT_cols)].shape[0]
        FP = only_PD.shape[0] # number of prediction sites without ground truth sites
        FN = only_GT.shape[0] # number of ground truth sites without prediction sites

        metric_name = f"Sensitivity:{window}nt"
        sensitivity = apa.sensitivity(TP, FN)
        metrics[metric_name] = [sensitivity, 0]

        metric_name = f"Precision:{window}nt"
        precision = apa.precision(TP, FP)
        metrics[metric_name] = [precision, 0]

        metric_name = f"F1_score:{window}nt"
        f1_score = apa.f1_score(precision, sensitivity)
        metrics[metric_name] = [f1_score, 0]

        metric_name = f"Jaccard_index:{window}nt"
        jaccard_index = apa.jaccard(TP, FP, FN)
        metrics[metric_name] = [jaccard_index, 0]

        ## Return-type dependent
        for return_df_type in all_return_df_types:
            if return_df_type == "all_GT":
                # add GT sites with no PD match
                sites = pd.concat([matched, only_GT])
            elif return_df_type == "union":
                # add GT sites with no PD match AND PD sites with no GT match
                sites = pd.concat([matched, only_GT, only_PD])
            elif return_df_type == "intersection":
                # do not consider any unmatched sites
                sites = matched
            else:
                raise ValueError(f"The variable return_df_type did not match any known string. Actual: {return_df_type}. Expected: all_GT, union or intersection.")
            
            sites.sort_values(by=['chrom_g', 'chromStart_g', 'chromEnd_g', 'chromStart_g'], inplace=True, ascending=[True, True, True, True])

            # METRIC: correlation coefficients
            #################################
            # Pearson correlation
            correlation_pearson = apa.corr_Pearson_with_gt(sites)
            # Key: exact name of metric as it appears in specification
            metric_name = f"Pearson_r:{return_df_type}:{window}nt"
            # Value: List of [variable_holding_metric, std_err]
            metrics[metric_name] = [correlation_pearson, 0]

            # Spearman correlation
            correlation_spearman = apa.corr_Spearman_with_gt(sites)
            # Key: exact name of metric as it appears in specification
            metric_name = f"Spearman_r:{return_df_type}:{window}nt"
            # Value: List of [variable_holding_metric, std_err]
            metrics[metric_name] = [correlation_spearman, 0]

            # METRIC: correlation coefficient of relative pas usage
            ####################
            # Only calculate metric for first window size.
            if window == windows[0]:
                normalised_df = apa.relative_pas_usage(sites, genome)

                # Pearson
                correlation_Pearson_rel_use = apa.corr_Pearson_with_gt(normalised_df)
                # Key: exact name of metric as it appears in specification
                metric_name = f"Pearson_r_relative:{return_df_type}:{window}nt"
                # Value: List of [variable_holding_metric, std_err]
                metrics[metric_name] = [correlation_Pearson_rel_use, 0]

                # Spearman
                correlation_Spearman_rel_use = apa.corr_Spearman_with_gt(normalised_df)
                # Key: exact name of metric as it appears in specification
                metric_name = f"Spearman_r_relative:{return_df_type}:{window}nt"
                # Value: List of [variable_holding_metric, std_err]
                metrics[metric_name] = [correlation_Spearman_rel_use, 0]

    # for the challenge, create all assessment json objects and append them to all_assessments
    for key, value in metrics.items():
        object_id = base_id + key
        assessment_object = JSON_templates.write_assessment_dataset(object_id, community, challenge, participant, key, value[0], value[1])
        all_assessments.append(assessment_object)

    # once all assessments have been added, print to json file
    with io.open(out_path,
                 mode='w', encoding="utf-8") as f:
        jdata = json.dumps(all_assessments, sort_keys=True, indent=4, separators=(',', ': '))
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
    parser.add_argument("-tpm", "--tpm_threshold", 
        help="Expression filter for predictions. PolyA sites with smaller or equal transcripts per million (tpm) will be removed before metric computation.", 
        required=True, type=int)

    args = parser.parse_args()


    main(args)
