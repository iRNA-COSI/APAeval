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
parser.add_argument("-c", "--challenge_ids", nargs='+', help="List of challenge ids selected by the user, separated by spaces", required=True)
parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
parser.add_argument("-o", "--output", help="output path where participant JSON file will be written", required=True)
parser.add_argument("-gtf", "--genome_dir", help="genome annotation directory. Used for relative PAS usage calculation. Directory needs to contain genome files with matching organism name from challenge.", required=True)

args = parser.parse_args()

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

def main(args):

    # input parameters
    participant_input = args.participant_data
    community = args.community_name
    challenges = args.challenge_ids
    participant_name = args.participant_name
    out_path = args.output
    genome_path = args.genome_dir
    gtf = select_genome_file(challenges[0], genome_path)
    gtf, with_chr = open(gtf,'r'), []
    for i in range(50):
        ln=gtf.readline()
        if not ln.startswith('#'):
            if ln.startswith('chr'):
                with_chr.append(True)
            else:
                with_chr.append(False)
    with_chr=list(set(with_chr)) ##should be only 1

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            print(os.path.dirname(out_path))
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a") : pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    validate_input_data(participant_input, community, challenges, participant_name, out_path, with_chr)



def  validate_input_data(participant_input, community, challenges, participant_name, out_path, with_chr):
    # get participant output (= input to be validated)
    try:
        participant_data = pandas.read_csv(participant_input, sep='\t',
                                           comment="#", header=0)
    except:
        sys.exit("ERROR: Submitted data file {} could not be read!".format(participant_input))
    
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
    if with_chr[0] == True:
        accepted_chr = ["chr"+str(i) for i in range(1,23)]
        accepted_chr.append('chrX')
        accepted_chr.append('chrY')
    else:
        accepted_chr = [str(i) for i in range(1,23)]
        accepted_chr.append('X')
        accepted_chr.append('Y')
    data_chr = list(set(participant_data.iloc[:, 0].values))
    chr_check = [str(chr) in accepted_chr for chr in data_chr].count(False) == 0
    
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
