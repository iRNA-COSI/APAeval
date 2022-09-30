#!/usr/bin/python3

import sys
import argparse
import csv

def parse_args(args=None):
    Description = "Reformat IsoSCM bed file into the output file of identification"
    Epilog = "Example usage: python postprcoess.py <FILE_IN> <RUN_IDENTIFICATION> <IDENTIFICATION_FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input IsoSCM output bed file.")
    parser.add_argument("RUN_MODE", help="Run mode could be 'identification'")
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
            # obtain proximal PAS
            chrom = row[2].split(":")[0]
            # need to do -1 to convert from 1-based gtf file coordinate to 0-based bed file coordinate
            start = int(row[2].split(":")[1].split("-")[1]) - 1
            end = start + 1
            strand = row[8]
            name = '|'.join([chrom, str(start) + ":" + str(end), strand])
            output = [chrom, str(start), str(end), name, ".", strand]
            outputs.append(output)

            # obtain distal PAS
            # if positive strand, distal PAS is the end of the downstream segment
            if strand == '+':
                start = int(row[6].split(":")[1].split("-")[1]) - 1
            # if negative strand, distal PAS is the end of the upstrea segment
            elif strand == '-':
                start = int(row[6]split(":")[1].split("-")[1]) - 1
            end = start + 1
            name = '|'.join([chrom, str(start) + ":" + str(end), strand])
            output = [chrom, str(start), str(end), name, ".", strand]
            outputs.append(output)

    return outputs


def write_identification_output(outputs, file_out):
    file_out = open(file_out, "wt")
    for output in outputs:
        file_out.write("\t".join(output) + "\n")
    file_out.close()


def main(args=None):
    args = parse_args(args)
    outputs = parse_input_file(args.FILE_IN)
    if args.RUN_MODE == 'identification':
        write_identification_output(outputs, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
