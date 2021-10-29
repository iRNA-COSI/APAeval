#!/usr/bin/env python3

import configparser
import argparse
import sys
import os

def main(base_dir,
         treat_dir,
         pas_base_dir,
         pas_treat_dir,
         annotation,
         genome,
         extended,
         all_events,
         output_dir,
         output_file):
    '''
    '''
    # Essentially a dict of dictionaries
    # Keys are the 'headers' (e.g. [INPUT_RNAseq])
    # Value = dictionary of parameter: value for each param in header/section
    config = configparser.ConfigParser()

    # Names of directories storing input files for each condition
    config["INPUT_RNAseq"] = {"input1": os.path.abspath(base_dir),
                              "input2": os.path.abspath(treat_dir)}

    # Names of directories storing PAS-seq input files for each condition
    if pas_base_dir is None:
        pas_base_dir = "NULL"

    else:
        pas_base_dir = os.path.abspath(pas_base_dir)

    if pas_treat_dir is None:
        pas_treat_dir = "NULL"

    else:
        pas_treat_dir = os.path.abspath(pas_treat_dir)


    config["INPUT_PASseq"] = {"pas1": pas_base_dir,
                              "pas2": pas_treat_dir}

    # Paths to GenePred annotation fle & genome fasta
    config["ANNOTATION"] = {"annotation": os.path.abspath(annotation),
                            "genome": os.path.abspath(genome)}


    config["Extended_3UTR"] = {"extended": extended}

    config["All_events"] = {"All": all_events}

    config["OUTPUT_FOLDER"] = {"output_dir": os.path.abspath(output_dir)}

    with open(output_file, "w") as configfile:
        config.write(configfile)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-b","--base-condition-dir",
                        required=True,
                        dest="base_cond_dir",
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to directory containing bam files for 'base' condition only")

    parser.add_argument("-t",
                        "--treatment-condition-dir",
                        required=True,
                        dest="treat_cond_dir",
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to directory containing bam files for 'treatment' condition only")

    parser.add_argument("--pas-seq-base-condition-dir",
                        dest="pas_seq_base_cond_dir",
                        default=None,
                        help="Path to directory containing PAS-seq BAM files for 'base' condition only")

    parser.add_argument("--pas-seq-treament-condition-dir",
                        dest="pas_seq_treat_cond_dir",
                        default=None,
                        help="Path to directory containing PAS-seq BAM files for 'treatment' condition only")

    parser.add_argument("-a",
                        "--annotation",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to input annotation file ('genePredExt' format)")

    parser.add_argument("-f",
                        "--fasta",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to input genome sequence fasta file")

    parser.add_argument("-e",
                        "--extended-mode",
                        choices=["yes","no"],
                        default="no",
                        type=str,
                        help="Whether to use APA-Scan's 'Extended-3'UTR' mode ('yes') or switch it off ('no')")

    parser.add_argument("--all-events",
                        choices=["yes","no"],
                        default="no",
                        type=str,
                        help="Whether to report all candidate cleavage sites of a gene whether significant or not. If 'no', most significant event for each gene will be reported")

    parser.add_argument("-d",
                        "--output-dir",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to output directory under which to store output files")

    parser.add_argument("-o", "--output-file",
                        required=True,
                        type=str,
                        default="configuration.ini",
                        help="Name of output APA-Scan configuration file")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.base_cond_dir,
         args.treat_cond_dir,
         args.pas_seq_base_cond_dir,
         args.pas_seq_treat_cond_dir,
         args.annotation,
         args.fasta,
         args.extended_mode,
         args.all_events,
         args.output_dir,
         args.output_file)
