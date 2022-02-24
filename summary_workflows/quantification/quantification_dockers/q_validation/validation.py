#!/usr/bin/env python3

from __future__ import division, print_function
import pandas
import numpy as np
import os, json
import sys
from argparse import ArgumentParser
from JSON_templates import JSON_templates

parser = ArgumentParser()
parser.add_argument("-i", "--participant_data", help="execution workflow prediction outputs", required=True)
parser.add_argument("-com", "--community_name", help="name of benchmarking community", required=True)
parser.add_argument("-c", "--challenge_types", nargs='+', help="list challenges selected by the user, separated by spaces", required=True)
parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
parser.add_argument("-o", "--output", help="output path where participant JSON file will be written",
                    required=True)

args = parser.parse_args()


def main(args):

    # input parameters
    input_participant = args.participant_data
    community = args.community_name
    challenges = args.challenge_types
    participant_name = args.participant_name
    out_path = args.output

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            print(os.path.dirname(out_path))
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a") : pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    validate_input_data(input_participant, community, challenges, participant_name, out_path)



def  validate_input_data(input_participant, community, challenges, participant_name, out_path):
    # get participant output (= input to be validated)
    try:
        participant_data = pandas.read_csv(input_participant, sep='\t',
                                           comment="#", header=0)
    except:
        sys.exit("ERROR: Submitted data file {} could not be read!".format(input_participant))
    
    #---------------------------------------------------
    # INPUT FILE VALIDATION
    # FOR APAeval QUANTIFICATION:
    # Check for valid bed6 format

    ## check number of columns
    n_col_check = len(participant_data.columns) == 6
    ## check start and end coordinates
    coord_check = participant_data.dtypes[1] == np.int64 and participant_data.dtypes[2] == np.int64
    ## check strands
    strands = list(set(participant_data.iloc[:, 5].values))
    strand_check = len(strands) == 2 and strands.count('-')+strands.count('+') == 2
    ## check ref seq format of chromosomes
    accepted_chr = [str(i) for i in range(1,23)]
    accepted_chr.append('X')
    accepted_chr.append('Y')
    data_chr = list(set(participant_data.iloc[:, 0].values))
    chr_check = [chr in accepted_chr for chr in data_chr].count(False) == 0
    
    ## All checks true?
    validated = False
    if n_col_check and coord_check and strand_check and chr_check:
        validated = True
    else:
        print("WARNING: Submitted data does not comply with required bed format.")
        validated = False
    #----------------------------------------------------


    data_id = community + ":" + participant_name
    output_json = JSON_templates.write_participant_dataset(data_id, community, challenges, participant_name, validated)

    # print file

    with open(out_path , 'w') as f:
        json.dump(output_json, f, sort_keys=True, indent=4, separators=(',', ': '))

    if validated == True:

        sys.exit(0)
    else:
        sys.exit("ERROR: Submitted data is not in valid format! Please check " + out_path)


if __name__ == '__main__':

    main(args)
