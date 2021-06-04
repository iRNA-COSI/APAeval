#!/usr/bin/env python3

import io
import json
import os
import fnmatch
from argparse import ArgumentParser


def main(args):
    # input parameters
    metrics_dir = args.metrics_data
    participant_dir = args.participant_data
    aggregation_dir = args.aggregation_data
    out_path = args.output

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    data_model_file = []

    data_model_file = join_json_files(participant_dir, data_model_file, "*.json")
    data_model_file = join_json_files(metrics_dir, data_model_file, "*.json")
    data_model_file = join_json_files(aggregation_dir, data_model_file, "*_summary.json")

    # write the merged data model file to json output
    with open(out_path, mode='w', encoding="utf-8") as f:
        json.dump(data_model_file, f, sort_keys=True, indent=4, separators=(',', ': '))

def join_json_files(data_directory, data_model_file, file_extension):


    # add minimal datasets to data model file
    if os.path.isfile(data_directory):
        with io.open(data_directory, mode='r', encoding="utf-8") as f:
            content = json.load(f)
            if isinstance(content, dict):
                data_model_file.append(content)
            else:
                data_model_file.extend(content)

    elif os.path.isdir(data_directory):  # if it is a directory loop over all files and search for all json

        for subdir, dirs, files in os.walk(data_directory):
            for file in files:
                abs_result_file = os.path.join(subdir, file)
                if fnmatch.fnmatch(abs_result_file, file_extension) and os.path.isfile(abs_result_file):
                    with io.open(abs_result_file, mode='r', encoding="utf-8") as f:
                        content = json.load(f)
                        if isinstance(content, dict):
                            data_model_file.append(content)
                        else:
                            data_model_file.extend(content)

    return data_model_file



if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-p", "--participant_data", help="path where the data for the participant is stored", required=True)
    parser.add_argument("-m", "--metrics_data", help="path where the data for the assessment metrics are stored", required=True)
    parser.add_argument("-a", "--aggregation_data", help="dir where the data for benchmark summary/aggregation are stored",
                        required=True)
    parser.add_argument("-o", "--output", help="output path where the minimal dataset JSON file will be written", required=True)

    args = parser.parse_args()

    main(args)

