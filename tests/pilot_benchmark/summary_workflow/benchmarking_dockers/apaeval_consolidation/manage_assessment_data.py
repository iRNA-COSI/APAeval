from __future__ import division
import io
import shutil
import json
import os
import fnmatch
from argparse import ArgumentParser
import numpy as np
from assessment_chart import assessment_chart

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

    with io.open(participant_path, mode='r', encoding="utf-8") as f:
        result = json.load(f)
        for item in result:
            participant_data.setdefault(item['challenge_id'], []).append (item)

    return participant_data

def generate_manifest(data_dir,output_dir, participant_data):

    info = []

    for cancer, metrics_file in participant_data.items():

        cancer_dir = os.path.join(output_dir,cancer)
        if not os.path.exists(cancer_dir):
            os.makedirs(cancer_dir)
        participants = []
        
        cancer_oeb_data = os.path.join(data_dir, cancer+".json")

        if os.path.isfile(cancer_oeb_data):
            # Transferring the public participants data
            with io.open(cancer_oeb_data, mode='r', encoding="utf-8") as f:
                aggregation_file = json.load(f)
                # get id for metrics in x and y axis
                metric_X = aggregation_file["datalink"]["inline_data"]["visualization"]["x_axis"]
                metric_Y = aggregation_file["datalink"]["inline_data"]["visualization"]["y_axis"]

                # add new participant data to aggregation file
                new_participant_data = {}
                for metrics_data in metrics_file:

                    participant_id = metrics_data["participant_id"]
                    if metrics_data["metrics"]["metric_id"] == metric_X:
                        new_participant_data ["metric_x"] = metrics_data["metrics"]["value"]
                    elif metrics_data["metrics"]["metric_id"] == metric_Y:
                        new_participant_data ["metric_y"] = metrics_data["metrics"]["value"]

                # copy the assessment files to output directory
                rel_new_location = participant_id + ".json"
                new_location = os.path.join(cancer_dir, rel_new_location)
                with open(new_location, 'w') as f:
                    json.dump(metrics_file, f, sort_keys=True, indent=4, separators=(',', ': '))

                new_participant_data["participant_id"] = participant_id
                aggregation_file["datalink"]["inline_data"]["challenge_participants"].append(new_participant_data)

                # add the rest of participants to manifest
                for name in aggregation_file["datalink"]["inline_data"]["challenge_participants"]:
                    participants.append(name["participant_id"])

                #copy the updated aggregation file to output directory
                summary_dir = os.path.join(cancer_dir,cancer + "_summary.json")
                with open(summary_dir, 'w') as f:
                    json.dump(aggregation_file, f, sort_keys=True, indent=4, separators=(',', ': '))



        # Let's draw the assessment charts!
        assessment_chart.print_chart(cancer_dir, summary_dir,cancer, "RAW")
        assessment_chart.print_chart(cancer_dir, summary_dir,cancer, "SQR")
        assessment_chart.print_chart(cancer_dir, summary_dir,cancer, "DIAG")

        #generate manifest
        obj = {
            "id" : cancer,
            "participants": participants
        }

        info.append(obj)

    with io.open(os.path.join(output_dir, "Manifest.json"), mode='w', encoding="utf-8") as f:
        jdata = json.dumps(info, f, sort_keys=True, indent=4, separators=(',', ': '))
        f.write(unicode(jdata,"utf-8"))


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-p", "--participant_data", help="path where the data for the participant is stored", required=True)
    parser.add_argument("-b", "--benchmark_data", help="dir where the data for the benchmark are stored", required=True)
    parser.add_argument("-o", "--output", help="output directory where the manifest and output JSON files will be written", required=True)

    args = parser.parse_args()

    main(args)
