#!/usr/bin/env python

import os
import sys
import argparse

def parse_args(args=None):
    Description = "Reformat APAeval samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context='Line', context_str=''):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != '' and context_str != '':
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(error, context.strip(), context_str.strip())
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq1,fastq2,bam,bai,gff,fasta,bed,mart_export
    """

    input_extensions = []
    sample_info_list = []
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 8
        HEADER = ['sample','fastq1','fastq2','bam', 'bai','gff','fasta','bed','mart_export']
        header = fin.readline().strip().split(",")
        if header[:len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip() for x in line.strip().split(",")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error("Invalid number of columns (minimum = {})!".format(len(HEADER)), 'Line', line)

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error("Invalid number of populated columns (minimum = {})!".format(MIN_COLS), 'Line', line)

            ## Check group name entries
            sample, fastq1, fastq2, bam, bai, gff, fasta, bed, mart_export = lspl[:len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Sample entry contains spaces!", 'Line', line)
            else:
                print_error("Sample entry has not been specified!", 'Line', line)

            ## Check fastq1 extension
            if fastq1:
                if fastq1.find(" ") != -1:
                    print_error("fastq1 contains spaces!", 'Line', line)
                if not fastq1.endswith(".fastq") and not fastq1.endswith(".fastq.gz") :
                    print_error("fastq1 does not have extension '.fastq' or '.fastq.gz'", 'Line', line)

            ## Check fastq2 extension
            if fastq2:
                if fastq2.find(" ") != -1:
                    print_error("fastq2 contains spaces!", 'Line', line)
                if not fastq2.endswith(".fastq") and not fastq2.endswith(".fastq.gz"):
                    print_error("fastq2 does not have extension '.fastq' or '.fastq.gz'", 'Line', line)

            ## Check bam extension
            if bam:
                if bam.find(" ") != -1:
                    print_error("bam contains spaces!", 'Line', line)
                if not bam.endswith(".bam"):
                    print_error("bam does not have extension 'bam'", 'Line', line)

            ## Check gff extension
            if gff:
                if gff.find(" ") != -1:
                    print_error("gff contains spaces!", 'Line', line)
                if not gff.endswith(".gff3"):
                    print_error("gff does not have extension '.gff3'", 'Line', line)

            ## Check fasta entries
            if fasta:
                if fasta.find(' ') != -1:
                    print_error("fasta entry contains spaces!",'Line', line)
                if len(fasta.split('.')) > 1:
                    if fasta[-6:] != '.fasta' and fasta[-3:] != '.fa' and fasta[-9:] != '.fasta.gz' and fasta[-6:] != '.fa.gz':
                        print_error("Genome entry does not have extension '.fasta', '.fa', '.fasta.gz' or '.fa.gz'!",'Line', line)

            ## Check 3UTR_bed extension
            if bed:
                if bed.find(" ") != -1:
                    print_error("bed contains spaces!", 'Line', line)
                if not bed.endswith(".bed"):
                    print_error("bed does not have extension '.bed'", 'Line', line)

            ## Create sample mapping dictionary = {group: {replicate : [ barcode, input_file, genome, gtf, is_transcripts ]}}
            sample_info = [ sample, fastq1, fastq2, bam, bai, gff, fasta, bed, mart_export ]
            sample_info_list.append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_info_list) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(['sample','fastq1','fastq2','bam', 'bai','gff','fasta','bed','mart_export']) + "\n")
            for sample_info in sample_info_list:
                ### Write to file
                fout.write(",".join(sample_info)+"\n")

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())

