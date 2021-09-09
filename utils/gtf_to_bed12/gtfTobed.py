#!/usr/bin/env python3

import os
# import sys
import argparse
import subprocess
# import glob

def main():
    parser = argparse.ArgumentParser(description='GTF to BED')
    parser.add_argument('--gtf', dest='gtf', type=str,
                        help='Input gtf file.', required=True)
    parser.add_argument('--out_bed', dest='out_bed', type=str,
                        help='Output BED12 file.', required=True)
    args = parser.parse_args()
    gtf=args.gtf
    bed=args.out_bed
    genePred=".".join(gtf.split(".")[:-1])+".genePred"
    #bed=".".join(gtf.split(".")[:-1])+".bed"
    cmd1 = ['/gtfToGenePred',  gtf, genePred]
    print(' '.join(cmd1))
    subprocess.run(' '.join(cmd1), shell=True, check=True)
    cmd2 = ['/genePredToBed', genePred, bed]
    print(' '.join(cmd2))
    subprocess.run(' '.join(cmd2), shell=True, check=True)
    os.remove(genePred)

if __name__ == "__main__":
    main()
