#!/usr/bin/env python3

from __future__ import division, print_function
import pandas
import numpy as np
import os, json
import sys
from argparse import ArgumentParser
import JSON_templates

parser = ArgumentParser()
parser.add_argument("-i", "--participant_data", help="Execution workflow output file to be validated", required=True)
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

    print(f"INFO: input {participant_input}")
    print(f"INFO: Possible challenges {challenges}")

    challenge = [c for c in challenges if c.split('.')[0] in str(participant_input)][0]
    
    print(f"INFO: Selected challenge {challenge}")

    # Get matching annotation file
    gtf = select_genome_file(challenge, genome_path)
    print(f"INFO: Selected genome file {gtf}")
    chr_names = list()

    with open(gtf, 'r') as f:
        for row in f:
            if not row.startswith('#'):
                # gtf is always tab-separated and first column is always seqname.
                seqname = row.split('\t')[0]
                if seqname not in chr_names:
                    chr_names.append(seqname)
    seqnames_wchr = [s for s in chr_names if 'chr' in s]
    assert len(seqnames_wchr) == len(chr_names) or len(seqnames_wchr) == 0, \
        f"WARNING: {genome_path} has a mix of chromosome name formats!"

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            print(os.path.dirname(out_path))
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a") : pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    validate_input_data(participant_input, community, challenge, participant_name, out_path, chr_names)



def  validate_input_data(infile, community, challenge, participant_name, out_path, chr_names):

    validated = False

    # get participant output (= input to be validated)
    try:
        participant_data = pandas.read_csv(infile, sep='\t',
                                        comment="#", header=None)
    except:
        sys.exit("ERROR: Submitted data file {} could not be read!".format(infile))
    
    #---------------------------------------------------
    # INPUT FILE VALIDATION
    # FOR APAeval RELATIVE QUANTIFICATION:
    # Check for valid bed6 format

    ## check number of columns
    n_col_check = len(participant_data.columns) == 6
    print(f"INFO: Columns check returned {n_col_check}")
    ## check start and end coordinates
    coord_check = participant_data.dtypes[1] == np.int64 and participant_data.dtypes[2] == np.int64
    print(f"INFO: Coordinate check returned {coord_check}")
    ## check strands
    strands = list(set(participant_data.iloc[:, 5].values))
    strand_check = len(strands) == 2 and strands.count('-')+strands.count('+') == 2
    print(f"INFO: Strand check returned {strand_check}")
    ## check ref seq format of chromosomes
    accepted_chr = chr_names
    data_chr = list(set(participant_data.iloc[:, 0].values))
    chr_check = [str(chr) in accepted_chr for chr in data_chr].count(False) == 0
    print(f"INFO: Chromosome check returned {chr_check}")
    
    ## All checks true?
    if n_col_check and coord_check and strand_check and chr_check:
        validated = True
    else:
        print(f"WARNING: Submitted file {infile} does not comply with required bed format.")
        validated = False
    #----------------------------------------------------


    data_id = community + ":" + participant_name
    output_json = JSON_templates.write_participant_dataset(data_id, community, challenge, participant_name, validated)

    # print validated participant file
    with open(out_path , 'w') as f:
        json.dump(output_json, f, sort_keys=True, indent=4, separators=(',', ': '))

    # Only pass if all input files are valid
    if validated == True:
        sys.exit(0)
    else:
        sys.exit("ERROR: One or more of the submitted files don't comply with APAeval specified format! Please check " + out_path)


if __name__ == '__main__':

    main(args)
