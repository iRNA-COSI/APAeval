
import csv
import time
import sys
from operator import itemgetter
import math
import pandas as pd
from Bio import SeqIO
import re
import bisect
from bisect import bisect_left
from scipy.stats import chisquare
import numpy as np
import methods
import glob, os
import preprocess
import configparser

speciesFlag, inputFlag, outFlag = 0, 0, 0

startTime = time.time()

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

config = configparser.ConfigParser()
config.read('configuration.ini')

input1_dir = config['INPUT_RNAseq']['input1']
input2_dir = config['INPUT_RNAseq']['input2']
if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]
pasSeq_dir1 = config['INPUT_PASseq']['pas1']
pasSeq_dir2 = config['INPUT_PASseq']['pas2']
if pasSeq_dir1[:-1] == "/":
	pasSeq_dir1 = pasSeq_dir1[:-1]
if pasSeq_dir2[:-1] == "/":
	pasSeq_dir2 = pasSeq_dir2[:-1]
output_dir = config['OUTPUT_FOLDER']['output_dir']
if output_dir[-1] != "/":
	output_dir += "/"
extended = config['Extended_3UTR']['extended']

os.makedirs(output_dir, exist_ok=True)
inp_annotation = config['ANNOTATION']['annotation']
ref_genome = config['ANNOTATION']['genome']
g1_name = input1_dir.split("/")[-1]
g2_name = input2_dir.split("/")[-1]
print(input1_dir)
print(input2_dir)
print("RNA-seq input 1 dir:", input1_dir)
print("RNA-seq input 2 dir:", input2_dir)
print("3'-end-seq input 1 dir:", pasSeq_dir1)
print("3'-end-seq input 2 dir:", pasSeq_dir2)
print("Output Dir:", output_dir)
print("Annotation:", inp_annotation, ref_genome, "\n\n")

print("Loading chromosomes...")
chromosomes = []
with open(inp_annotation, 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter='\t')
    headers = next(f)
    annotList = list(reader)
    for rows in annotList:
    	if '_' not in rows[2] and rows[2]!='chrM':
	    	chromosomes.append(rows[2])
    chr_set = set(chromosomes)
    chromosomes = list(chr_set)

chromDict = methods.makeChromDict(chromosomes, inp_annotation)

print("Creating read coverage files for RNA-seq data...")
os.chdir(input1_dir)
for sample1 in glob.glob("*.bam"):
    preprocess.SamtoText(input1_dir, sample1, chromosomes)
os.chdir(input2_dir)
for sample2 in glob.glob("*.bam"):
    preprocess.SamtoText(input2_dir, sample2, chromosomes)

result_filename = output_dir+"APA_Scan_"+g1_name+"_Vs_"+g2_name
if pasSeq_dir1 == 'NULL' or pasSeq_dir2=='NULL':
	s1_namelist = list_dirs(input1_dir)
	s2_namelist = list_dirs(input2_dir)

	print("Preparing result using RNA-seq data only")
	methods.Get_Signal_Positions(chromosomes, chromDict, inp_annotation, ref_genome, output_dir)  # EDIT HERE
    # DROPPED THE LAST (6TH) ARGUMENT "extnede" FROM THIS FUNCTION AS I GOT:
    #   "TypeError: Get_Signal_Positions() takes 5 positional arguments but 6 were given"
	methods.with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, result_filename)

else:
	print("Creating read coverage files for 3'-end-seq data...")
	os.chdir(pasSeq_dir1)
	for sample1 in glob.glob("*.bam"):
	    preprocess.SamtoText(pasSeq_dir1, sample1, chromosomes)
	os.chdir(pasSeq_dir2)
	for sample2 in glob.glob("*.bam"):
	    preprocess.SamtoText(pasSeq_dir2, sample2, chromosomes)

	p1_namelist = list_dirs(pasSeq_dir1)
	p2_namelist = list_dirs(pasSeq_dir2)
	p1_name = pasSeq_dir1.split("/")[-1]
	p2_name = pasSeq_dir2.split("/")[-1]

	filename = output_dir+'PA_peak_positions.csv'
	#methods.Get_Peak_Positions(filename, chromosomes, inp_annotation, pasSeq_dir1, pasSeq_dir2, p1_name, p2_name, output_dir)
	#methods.with_PA_peaks(chromosomes, input1_dir, input2_dir, g1_name, g2_name, filename, output_dir, result_filename)
	methods.Generate_withPasSeqData(filename, chromosomes, inp_annotation, pasSeq_dir1, pasSeq_dir2, p1_name, p2_name, output_dir)
	methods.Quantification(chromosomes, input1_dir, input2_dir, g1_name, g2_name, filename, result_filename)

print("Total time:", round((time.time() - startTime)/60, 2), "minutes.")
read_file = pd.read_csv(result_filename+".csv", delimiter = '\t')
read_file.to_excel (result_filename+".xlsx", index = None, header=True)
os.remove(result_filename+".csv")
