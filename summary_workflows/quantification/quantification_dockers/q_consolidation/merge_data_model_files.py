#!/usr/bin/env python3

import io
import json
import os
import fnmatch
from argparse import ArgumentParser


def main(args):
    # input parameters
    metrics_data = args.metrics_data
    validation_data = args.validation_data
    aggregation_data = args.aggregation_data
    challenge = args.challenge_id
    out_path = args.output

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)
            
    # This is the final consolidated output
    data_model_file = []

    # get the output files from previous steps and concatenate
    # from validation ("validated_participant_data")
    data_model_file = join_json_files(validation_data, data_model_file, "*.json")
    # from metrics ("assessment_out")
    data_model_file = join_json_files(metrics_data, data_model_file, "*.json")
    # from consolidation part 1 (manage_assessment_data.py), "sample_out/results/challenge.json"
    data_model_file = join_json_files(aggregation_data, data_model_file, "*" + challenge + ".json")

    # write the merged data model file to json output
    with open(out_path, mode='w', encoding="utf-8") as f:
        json.dump(data_model_file, f, sort_keys=True, indent=4, separators=(',', ': '))

def join_json_files(data_directory, data_model_file, file_extension):
    '''Add contents of specified file(s) to given json file
    Input:
    data_directory: Directory containing the file(s) to be added
    data_model_file: list to be extended with json objects (dicts)
    file_extension: filename or pattern to be matched, of the file(s) to be added
    Returns:
    data_model_file: data_model_file (which is a list) from input, extended by the json objects from the input file(s)
    '''

    # add minimal datasets to data model file
    if os.path.isfile(data_directory):
        with io.open(data_directory, mode='r', encoding="utf-8") as f:
            content = json.load(f)
            if isinstance(content, dict):
                data_model_file.append(content)
            else:
                data_model_file.extend(content)

    elif os.path.isdir(data_directory):  # if it is a directory loop over all files and search for the file with the given "extension" (=pattern, can be name)

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
    parser.add_argument("-v", "--validation_data", help="path to validated_participant_data.json", required=True)
    parser.add_argument("-m", "--metrics_data", help="path to assessment_datasets.json", required=True)
    parser.add_argument("-a", "--aggregation_data", help="dir where the data for benchmark summary/aggregation are stored",
                        required=True)
    parser.add_argument("-c", "--challenge_id", help="Id of the challenge", required=True)
    parser.add_argument("-o", "--output", help="output path where the minimal dataset JSON file will be written", required=True)

    args = parser.parse_args()

    main(args)

