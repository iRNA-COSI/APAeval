#!/usr/bin/python3

import sys
import argparse
import csv

def parse_args(args=None):
    Description = "Reformat IsoSCM bed file into the output file of identification and relative usage quantification challenge"
    Epilog = "Example usage: python postprcoess_identification.py <FILE_IN> <RUN_IDENTIFICATION> <RUN_RELATIVE_USAGE_QUANTIFICATION> <IDENTIFICATION_FILE_OUT> <RELATIVE_USAGE_QUANTIFICATION_FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input IsoSCM output bed file.")
    parser.add_argument("RUN_MODE", help="Run mode could be either 'identification' or 'relative_usage_quantification'")
    parser.add_argument("FILE_OUT", help="Name of output file")
    return parser.parse_args(args)


def parse_input_file(file_in):
    """
    This function parses IsoSCM output file to file for
    identification or relative usage quantification challenge
    :param file_in: IsoSCM compare step output file to be reformatted to relative usage quantification bed file
    :return: N/A
    """
    outputs = []
    with open(file_in) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        next(fd, None)
        for row in rd:
            chrom = row[2].split(":")[0]
            # need to do -1 to convert from 1-based gtf file coordinate to 0-based bed file coordinate
            start = int(row[2].split(":")[1].split("-")[1]) - 1
            end = start + 1
            strand = row[8]
            name = '|'.join([chrom, str(start) + ":" + str(end), strand])
            # score is taken from site usage column that has the site usage for both samples
            # since the two samples for IsoSCM compare are the same, the site usage values
            # are the same (e.g. 0.9895184913950934,0.9895184913950934), so we can arbitrarily take the first 
            score = row[11].split(",")[0]
            output = [chrom, str(start), str(end), name, score, strand]
            outputs.append(output)

    return outputs


def write_identification_output(outputs, file_out):
    file_out = open(file_out, "wt")
    for output in outputs:
        output[4] = "."
        file_out.write("\t".join(output) + "\n")
    file_out.close()


def write_relative_usage_quantification_output(outputs, file_out):
    file_out = open(file_out, "wt")
    for output in outputs:
        file_out.write("\t".join(output) + "\n")
    file_out.close()


def main(args=None):
    args = parse_args(args)
    outputs = parse_input_file(args.FILE_IN)
    if args.RUN_MODE == 'identification':
        write_identification_output(outputs, args.FILE_OUT)
    if args.RUN_MODE == 'relative_usage_quantification':
        write_relative_usage_quantification_output(outputs, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
