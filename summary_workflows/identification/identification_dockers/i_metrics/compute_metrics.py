#!/usr/bin/env python3

from __future__ import division
import io
import os
import json
from argparse import ArgumentParser
from JSON_templates import JSON_templates
from identification_metrics import identification_metrics as im

def main(args):

    # input parameters
    participant_input = args.participant_data
    gold_standards_dir = args.gold_standards_dir
    challenge_ids = args.challenge_ids
    participant = args.participant_name
    community = args.community_name
    out_path = args.output
    windows = args.windows

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    compute_metrics(participant_input, gold_standards_dir, challenge_ids, participant, community, out_path, windows)


def compute_metrics(participant_input, gold_standards_dir, challenge_ids, participant, community, out_path, windows):

    # define array that will hold the full set of assessment datasets
    all_assessments = []

    for challenge in challenge_ids:
        
        # ID prefix for assessment objects
        base_id = f"{community}:{challenge}_{participant}_"
        # Dict to store metric names and corresponding variables + stderr (which is currently not computed and set to 0)
        metrics = {}
        # ground truth file
        gold_standard = os.path.join(gold_standards_dir, challenge + ".bed")

        for window in windows:

            # METRIC: Precision
            ########################

            precision = im.precision(participant_input,gold_standard,window)

            # Key: exact name of metric as it appears in specification
            metric_name = f"Precision_{window}nt"
            # Value: List of [variable_holding_metric, std_err]
            metrics[metric_name] = [precision, 0]


            # METRIC: Sensitivity
            #################################

            sensitivity = im.sensitivity(participant_input,gold_standard,window)

            # Key: exact name of metric as it appears in specification
            metric_name = f"Sensitivity_{window}nt"
            # Value: List of [variable_holding_metric, std_err]
            metrics[metric_name] = [sensitivity, 0]
                    

            # METRIC: Multi-matched
            ####################

            multi = im.multi_matched(participant_input,gold_standard,window)
            
            # Key: exact name of metric as it appears in specification
            metric_name = f"Multi-matched_{window}nt"
            # Value: List of [variable_holding_metric, std_err]
            metrics[metric_name] = [multi, 0]

        # Outside of window loop:

        # METRIC: AUC
        ######################

        auc = im.auc(metrics)
            
        # Key: exact name of metric as it appears in specification
        metric_name = f"AUC"
        # Value: List of [variable_holding_metric, std_err]
        metrics[metric_name] = [auc, 0]

        # METRIC: Proportion in 3UTR
        #############################

        prop_utr = im.proportion_3UTR(participant_input, gtf)
            
        # Key: exact name of metric as it appears in specification
        metric_name = f"Proportion_in_3UTR"
        # Value: List of [variable_holding_metric, std_err]
        metrics[metric_name] = [prop_utr, 0]

        # METRIC: Genes correct PAS
        ################################

        correct_PAS = im.genes_correct_PAS(participant_input, gtf)
            
        # Key: exact name of metric as it appears in specification
        metric_name = f"Genes_correct_PAS"
        # Value: List of [variable_holding_metric, std_err]
        metrics[metric_name] = [correct_PAS, 0]

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


if __name__ == '__main__':
    
    parser = ArgumentParser()
    parser.add_argument("-i", "--participant_data", help="execution workflow prediction outputs", required=True)
    parser.add_argument("-c", "--challenge_ids", nargs='+', help="List of challenge ids selected by the user, separated by spaces", required=True)
    parser.add_argument("-g", "--gold_standards_dir", help="dir that contains gold standard datasets for current challenge", required=True)
    parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
    parser.add_argument("-com", "--community_name", help="name/id of benchmarking community", required=True)
    parser.add_argument("-o", "--output", help="output path where assessment JSON files will be written", required=True)
    parser.add_argument("-w", "--windows", nargs='+', help="windows (nt) for scanning for poly(A) sites; several window sizes separated by spaces", required=True, type=int)
    
    args = parser.parse_args()

    
    main(args)




