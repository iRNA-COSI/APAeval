#!/usr/bin/env python3

from __future__ import division
import io
import os
import pandas
import math
import json
from argparse import ArgumentParser
from JSON_templates import JSON_templates
from Leo_MatchPAsites_GT import Leo_MatchPAsites_GT

def main(args):

    # input parameters
    input_participant = args.participant_data
    gold_standards_dir = args.metrics_ref
    challenge_types = args.challenge_types
    participant = args.participant_name
    community = args.community_name
    out_path = args.output

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    compute_metrics(input_participant, gold_standards_dir, challenge_types, participant, community, out_path)

def convertRunTimeToSec(runtime_str):
    runtime_broken_down = runtime_str.split()
    if len(runtime_broken_down) == 1:
        runtime = int(runtime_broken_down[-1].strip('s'))
    elif len(runtime_broken_down) == 2:
        runtime = int(runtime_broken_down[-1].strip('s'))+int(runtime_broken_down[-2].strip('m'))*60
    elif len(runtime_broken_down) == 3:
        runtime = int(runtime_broken_down[-1].strip('s'))+int(runtime_broken_down[-2].strip('m'))*60+int(runtime_broken_down[-3].strip('h'))*60*60
    return runtime

def compute_metrics(input_participant, gold_standards_dir, challenge_types, participant, community, out_path):

    # get participant dataset
#    participant_data = pandas.read_csv(input_participant, sep='\t',
#                                       comment="#", header=0, index_col=0)

    # define array that will hold the full set of assessment datasets
    ALL_ASSESSMENTS = []

    for challenge in challenge_types:
        # get metrics dataset
#        metrics_data = pandas.read_csv(os.path.join(gold_standards_dir, challenge + ".txt"),
#                                       comment="#", header=0, sep='\t',index_col=0)
        metrics_data=os.path.join(gold_standards_dir, "siControl_R2.MACEseq.mm10.bed") ##PLEASE CHANGE THIS LATER

        # metric on the number of matched sites
        window = 15
        match_with_gt_run = Leo_MatchPAsites_GT.match_wtih_gt(input_participant,metrics_data,window)
        merged_bed_df, n_matched_sites, n_unmatched_sites = match_with_gt_run[0], match_with_gt_run[1], match_with_gt_run[2]

        # metric on correlation coffecient
        correlation=Leo_MatchPAsites_GT.corr_with_gt(merged_bed_df)

        #assessment_data = {'toolname': participant, 'x': TPR, 'y': acc, 'e': 0, 'challenge_type': challenge} #not used anywhere

        # get json assessment file for both metrics
        data_id_1 = community + ":" + challenge + "_runtime_" + participant + "_A"
        std_error= 0
        assessment_matched_sites = JSON_templates.write_assessment_dataset(data_id_1, community, challenge, participant, "n_matched_sites", n_matched_sites, std_error)

        data_id_2 = community + ":" + challenge + "_runtime_" + participant + "_A"
        std_error= 0
        assessment_unmatched_sites = JSON_templates.write_assessment_dataset(data_id_2, community, challenge, participant, "n_unmatched_sites", n_unmatched_sites, std_error)

        data_id_3 = community + ":" + challenge + "_memory_" + participant + "_A"
        std_error= 0
        assessment_correlation = JSON_templates.write_assessment_dataset(data_id_3, community, challenge, participant, "correlation", correlation, std_error)

        # push the two assessment datasets to the main dataset array
        ALL_ASSESSMENTS.extend([assessment_matched_sites, assessment_unmatched_sites, assessment_correlation])

    # once all assessments have been added, print to json file
    with io.open(out_path,
                 mode='w', encoding="utf-8") as f:
        jdata = json.dumps(ALL_ASSESSMENTS, sort_keys=True, indent=4, separators=(',', ': '))
        f.write(jdata)


if __name__ == '__main__':
    
    parser = ArgumentParser()
    parser.add_argument("-i", "--participant_data", help="list of challenge genes prediction", required=True)
    parser.add_argument("-c", "--challenge_types", nargs='+', help="list of types of challenge selected by the user, separated by spaces", required=True)
    parser.add_argument("-m", "--metrics_ref", help="dir that contains metrics reference datasets for all challenge types", required=True)
    parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
    parser.add_argument("-com", "--community_name", help="name/id of benchmarking community", required=True)
    parser.add_argument("-o", "--output", help="output path where assessment JSON files will be written", required=True)
    
    args = parser.parse_args()

    
    main(args)




