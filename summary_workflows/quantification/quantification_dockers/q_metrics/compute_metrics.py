#!/usr/bin/env python3

from __future__ import division
import io
import os
import json
import pandas as pd
from argparse import ArgumentParser
from JSON_templates import JSON_templates
from matchPAS import matchPAS

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

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    compute_metrics(participant_input, gold_standards_dir, challenge_ids, participant, community, out_path, windows, genome_path)


def compute_metrics(participant_input, gold_standards_dir, challenge_ids, participant, community, out_path, windows, genome_path):

    # define array that will hold the full set of assessment datasets
    all_assessments = []
    # define list with keywords for return type of dataframe from match_with_gt()
    all_return_df_types = ["all_GT", "union", "intersection"]

    for challenge in challenge_ids:

        # ID prefix for assessment objects
        base_id = f"{community}:{challenge}:{participant}:"
        # Dict to store metric names and corresponding variables + stderr (which is currently not computed and set to 0)
        metrics = {}
        # ground truth file
        gold_standard = os.path.join(gold_standards_dir, challenge + ".bed")
        # genome annotation file
        genome_file = select_genome_file(challenge, genome_path)
        # log to stdout
        print(f"INFO: In challenge {challenge}. Using genome file: {genome_file}")
        genome = matchPAS.load_genome(genome_file)

        for window in windows:
            # obtain dataframes with matched and unmatched sites
            matched, only_PD, only_GT = matchPAS.match_with_gt(f_PD=participant_input, f_GT=gold_standard, window=window)

            # METRIC: Expression unmatched sites
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
            
           
            #??METRIC: Percent expression unmatched
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
            # number of ground truth sites with matching prediction site and TPM > 0
            TP = matched.loc[matched["score_p"] > 0].shape[0]
            FP = only_PD.shape[0] # number of prediction sites without ground truth sites
            FN = only_GT.shape[0] # number of ground truth sites without prediction sites

            metric_name = f"Sensitivity:{window}nt"
            sensitivity = TP / (TP + FN)
            metrics[metric_name] = [sensitivity, 0]

            metric_name = f"Precision:{window}nt"
            precision = TP / (TP + FP)
            metrics[metric_name] = [precision, 0]

            metric_name = f"F1_score:{window}nt"
            f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)
            metrics[metric_name] = [f1_score, 0]

            metric_name = f"Jaccard_index:{window}nt"
            jaccard_index = TP / (TP + FP + FN)
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
                correlation_pearson = matchPAS.corr_Pearson_with_gt(sites)
                # Key: exact name of metric as it appears in specification
                metric_name = f"Pearson_r:{return_df_type}:{window}nt"
                # Value: List of [variable_holding_metric, std_err]
                metrics[metric_name] = [correlation_pearson, 0]

                # Spearman correlation
                correlation_spearman = matchPAS.corr_Spearman_with_gt(sites)
                # Key: exact name of metric as it appears in specification
                metric_name = f"Spearman_r:{return_df_type}:{window}nt"
                # Value: List of [variable_holding_metric, std_err]
                metrics[metric_name] = [correlation_spearman, 0]

                # METRIC: correlation coefficient of relative pas usage
                ####################
                # Only calculate metric for first window size.
                if window == windows[0]:
                    normalised_df = matchPAS.relative_pas_usage(sites, genome)

                    # Pearson
                    correlation_Pearson_rel_use = matchPAS.corr_Pearson_with_gt(normalised_df)
                    # Key: exact name of metric as it appears in specification
                    metric_name = f"Pearson_r_relative:{return_df_type}:{window}nt"
                    # Value: List of [variable_holding_metric, std_err]
                    metrics[metric_name] = [correlation_Pearson_rel_use, 0]

                    # Spearman
                    correlation_Spearman_rel_use = matchPAS.corr_Spearman_with_gt(normalised_df)
                    # Key: exact name of metric as it appears in specification
                    metric_name = f"Spearman_r_relative:{return_df_type}:{window}nt"
                    # Value: List of [variable_holding_metric, std_err]
                    metrics[metric_name] = [correlation_Spearman_rel_use, 0]

        # for each challenge, create all assessment json objects and append them to all_assessments
        for key, value in metrics.items():
            object_id = base_id + key
            assessment_object = JSON_templates.write_assessment_dataset(object_id, community, challenge, participant, key, value[0], value[1])
            all_assessments.append(assessment_object)

    # once all assessments have been added, print to json file
    with io.open(out_path,
                 mode='w', encoding="utf-8") as f:
        jdata = json.dumps(all_assessments, sort_keys=True, indent=4, separators=(',', ': '))
        f.write(jdata)


def select_genome_file(file_name, genome_path):
    """Select the genome file according to the organism.

    Requires that the file_name contains an expression containing organism
    information, which will be matched against the genome_path directory.
    The format should be: name.mm10.ext or name.hg38extension.ext, with
    matching genome annotations: gencode.mm10.gtf and gencode.hg38extension.gtf.
    Note: no check for the extension (e.g. gtf) is done.

    Args:
        file_name (str): Name containing organism information. Supported: mm* and hg*.
        genome_path (str): directory containing genome annotations in gtf format.

    Returns:
        str with genome file path.
    """
    GENOME_STRINGS = ["mm", "hg"]
    SPLITSTRING = "."
    assert os.path.exists(genome_path), f"Genome annotation directory not found: {genome_path}"
    file_components =  file_name.split(SPLITSTRING)
    # search for genome
    for genome_string in GENOME_STRINGS:
        match = [comp for comp in file_components if genome_string in comp]
        if len(match) != 0:
            break
    if len(match) == 0:
        raise ValueError(f"No genome string: {GENOME_STRINGS} in file_name: {file_name} found.")
    # find all genome files in genome_path
    for f in os.listdir(genome_path):
        # find exact match in file
        genome_match = [f for comp in f.split(SPLITSTRING) if match[0] == comp]
        if len(genome_match) != 0:
            break
    if len(genome_match) == 0:
        raise ValueError(f"No genome string: {GENOME_STRINGS} in genome_path: {genome_path} found.")
    # return file
    return os.path.join(genome_path, genome_match[0])


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-i", "--participant_data", help="execution workflow prediction outputs", required=True)
    parser.add_argument("-c", "--challenge_ids", nargs='+', help="List of challenge ids selected by the user, separated by spaces", required=True)
    parser.add_argument("-g", "--gold_standards_dir", help="dir that contains gold standard datasets for current challenge", required=True)
    parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
    parser.add_argument("-com", "--community_name", help="name/id of benchmarking community", required=True)
    parser.add_argument("-o", "--output", help="output path where assessment JSON files will be written", required=True)
    parser.add_argument("-w", "--windows", nargs='+', help="windows (nt) for scanning for poly(A) sites; several window sizes separated by spaces. The first entry is used for MSE calculation.", required=True, type=int)
    parser.add_argument("-gtf", "--genome_dir", help="genome annotation directory. Used for relative PAS usage calculation. Directory needs to contain genome files with matching organism name from challenge.", required=True)

    args = parser.parse_args()


    main(args)
