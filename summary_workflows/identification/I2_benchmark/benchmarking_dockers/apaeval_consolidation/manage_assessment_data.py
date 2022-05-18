#!/usr/bin/env python3

from __future__ import division
import shutil
import json
import os
import fnmatch
from argparse import ArgumentParser
import numpy as np
from assessment_chart import assessment_chart

DEFAULT_eventMark = '2018-04-05'

def main(args):

    # input parameters
    data_dir = args.benchmark_data
    participant_path = args.participant_data
    output_dir = args.output
    
    # Assuring the output directory does exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # read participant metrics
    participant_data = read_participant_data(participant_path)
    generate_manifest(data_dir, output_dir, participant_data)


def read_participant_data(participant_path):
    participant_data = {}

    with open(participant_path, mode='r', encoding="utf-8") as f:
        result = json.load(f)
        for item in result:
            participant_data.setdefault(item['challenge_id'], []).append (item)

    return participant_data

def generate_manifest(data_dir,output_dir, participant_data):

    info = []

    for challenge, metrics_file in participant_data.items():

        challenge_dir = os.path.join(output_dir,challenge)
        if not os.path.exists(challenge_dir):
            os.makedirs(challenge_dir)
        participants = []
        
        
        challenge_oeb_data_dir = os.path.join(data_dir, challenge)
        challenge_oeb_data = challenge_oeb_data_dir + ".json"

        if os.path.isfile(challenge_oeb_data):
            # Transferring the public participants data
            with open(challenge_oeb_data, mode='r', encoding="utf-8") as f:
                aggregation_file = json.load(f)
            
            # get default id for metrics in x and y axis
            metric_X = aggregation_file["datalink"]["inline_data"]["visualization"]["x_axis"]
            metric_Y = aggregation_file["datalink"]["inline_data"]["visualization"]["y_axis"]
        else:
            challenge_participants = []
            
            # get default id for metrics in x and y axis
            metric_X = None
            metric_Y = None
            for metrics_data in metrics_file:
                if metric_X is None:
                    metric_X = metrics_data["metrics"]["metric_id"]
                elif metric_Y is None:
                    metric_Y = metrics_data["metrics"]["metric_id"]
                else:
                    break
            
            # Setting the defaults in case nothing was found
            if metric_X is None:
                metric_X = "TPR"
            if metric_Y is None:
                metric_Y = "precision"
            
            aggregation_file = {
                "_id": "TCGA:{}_{}_Aggregation".format(DEFAULT_eventMark, challenge),
                "challenge_ids": [
                    challenge
                ],
                "datalink": {
                    "inline_data": {
                        "challenge_participants": challenge_participants,
                        "visualization": {
                            "type": "2D-plot",
                            "x_axis": metric_X,
                            "y_axis": metric_Y
                        }
                    }
                },
                "type": "aggregation"
            }
            
            # Get the info from the files in the directory
            if os.path.isdir(challenge_oeb_data_dir):
                print("Reading {}".format(challenge_oeb_data_dir))
                for entry in os.scandir(challenge_oeb_data_dir):
                    if entry.is_file() and entry.name.endswith(".json"):
                        with open(entry.path, mode="r", encoding="utf-8") as ep:
                            metrics_content = json.load(ep)
                            if metrics_content.get("challenge_type") == challenge:
                                challenge_participants.append({
                                    "metric_x": metrics_content["x"],
                                    "metric_y": metrics_content["y"],
                                    "participant_id": metrics_content["toolname"]
                                })

        # add new participant data to aggregation file
        new_participant_data = {}
        participant_id = '(unknown)'
        for metrics_data in metrics_file:
            participant_id = metrics_data["participant_id"]
            if metrics_data["metrics"]["metric_id"] == metric_X:
                new_participant_data["metric_x"] = metrics_data["metrics"]["value"]
            elif metrics_data["metrics"]["metric_id"] == metric_Y:
                new_participant_data["metric_y"] = metrics_data["metrics"]["value"]

        # copy the assessment files to output directory
        rel_new_location = participant_id + ".json"
        new_location = os.path.join(challenge_dir, rel_new_location)
        with open(new_location, mode='w', encoding="utf-8") as f:
            json.dump(metrics_file, f, sort_keys=True, indent=4, separators=(',', ': '))

        new_participant_data["participant_id"] = participant_id
        aggregation_file["datalink"]["inline_data"]["challenge_participants"].append(new_participant_data)

        # add the rest of participants to manifest
        for name in aggregation_file["datalink"]["inline_data"]["challenge_participants"]:
            participants.append(name["participant_id"])

        #copy the updated aggregation file to output directory
        summary_dir = os.path.join(challenge_dir,challenge + "_summary.json")
        with open(summary_dir, 'w') as f:
            json.dump(aggregation_file, f, sort_keys=True, indent=4, separators=(',', ': '))

        # Let's draw the assessment charts!
        assessment_chart.print_chart(challenge_dir, summary_dir, challenge, "RAW")
        assessment_chart.print_chart(challenge_dir, summary_dir, challenge, "SQR")
        assessment_chart.print_chart(challenge_dir, summary_dir, challenge, "DIAG")

        #generate manifest
        obj = {
            "id" : challenge,
            "participants": participants
        }

        info.append(obj)

    with open(os.path.join(output_dir, "Manifest.json"), mode='w', encoding="utf-8") as f:
        json.dump(info, f, sort_keys=True, indent=4, separators=(',', ': '))


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-p", "--participant_data", help="path where the data for the participant is stored", required=True)
    parser.add_argument("-b", "--benchmark_data", help="dir where the data for the benchmark are stored", required=True)
    parser.add_argument("-o", "--output", help="output directory where the manifest and output JSON files will be written", required=True)

    args = parser.parse_args()

    main(args)
