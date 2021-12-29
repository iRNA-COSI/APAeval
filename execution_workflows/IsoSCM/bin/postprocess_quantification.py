'''
cat ${ASSEMDIR}/${SAMPLE}.coverage.gtf | awk '$9 ~ /coverage/' | \
    sed -E 's/(.*coverage ")([^"]*)(.*)/\1\2\3\t\2/g' > tmp.gtf

awk -F"\t" -vOFS="\t" 'NR==FNR{bed[$1"_"$2]=$0;next} { for (b in bed) \
    {split(b,pos,"_"); if (pos[1]==$1 && pos[2]>=$4 && pos[2]<=$5) \
    {print bed[b]"\t"$10; delete bed[b]}} }' ${SAMPLE}_IsoSCM_01.bed tmp.gtf | \
    awk -vOFS="\t" '{print $1,$2,$3,$4,$7,$6}' > ${SAMPLE}_IsoSCM_02.bed
'''
#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

def parse_args(args=None):
    Description = "Reformat IsoSCM bed file into the output files of quantification and quantification challenges"
    Epilog = "Example usage: python postprcoess_quantification.py <FILE_IN> <QUANTIFICATION_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input IsoSCM output bed file.")
    parser.add_argument("QUANTIFICATION_OUT", help="Name of output file for quantification challenge")
    return parser.parse_args(args)


def reformat_bed(file_in, quantification_out):
    """
    This function reformats IsoSCM output file to file for
    quantification challenge
    :param file_in: txt file to be reformatted
    :return: N/A
    """
    quantification_out = open(quantification_out, "wt")

    df = pd.read_csv(file_in, sep='\t')
    for index, row in df.iterrows():
        chrom = row[0]
        start = row[3]
        end = int(start) + 1
        strand = row[6]
        name = '|'.join([chrom, str(start)+":"+str(end), strand])
        score = "."
        output = (chrom, start, end, name, score, strand)
        quantification_out.write("\t".join(output) + "\n")

    quantification_out.close()


def main(args=None):
    args = parse_args(args)
    reformat_bed(args.FILE_IN, args.QUANTIFICATION_OUT)


if __name__ == '__main__':
    sys.exit(main())
