#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import glob

def main():
    parser = argparse.ArgumentParser(description='GTF to BED')
    parser.add_argument('--gtf', dest='gtf', type=str,
                        help='Input gtf file.', required=True)
    args = parser.parse_args()
    gtf=args.gtf
    genePred=".".join(gtf.split(".")[:-1])+".genePred"
    bed=".".join(gtf.split(".")[:-1])+".bed"
    cmd1 = ['/app/gtfToGenePred',  gtf, genePred]
    print(' '.join(cmd1))
    subprocess.run(' '.join(cmd1), shell=True, check=True)
    cmd2 = ['/app/genePredToBed', genePred, bed]
    print(' '.join(cmd2))
    subprocess.run(' '.join(cmd2), shell=True, check=True)

if __name__ == "__main__":
    main()
