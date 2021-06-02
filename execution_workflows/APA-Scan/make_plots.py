import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import bisect
from bisect import bisect_left
import csv
import methods
import time
import sys
import os, glob
import configparser
import xlsxwriter
from matplotlib import rcParams
plt.rc('legend',**{'fontsize':10})

rcParams.update({
    'font.family':'arial',
    })
y_limit, y_limit2 = 0, 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

def Generate_read_coverate_plot(ax, pathin, sample, labelname, chrom, geneID, startAll, endAll, n1, p1):
	bam_file_reader= open(pathin+'/'+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	pos1 = bi_contains(position_row, startAll)
	pos2 = bi_contains(position_row, endAll)
	if(int(bam_list[pos2][1]) != endAll):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = endAll - startAll + 1
	
	for t in range(length):
		p.append(t+startAll)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][1])
		read = int(bam_list[t][2])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	if p1 == 0:
		labelname = labelname+" (RNA-seq)"
		global y_limit
		m = max(c)
		if m > y_limit:
			y_limit = m
		yy = y_limit
	else:
		labelname = labelname+" (3'-end-seq)"
		global y_limit2
		m = max(c)
		if m > y_limit2:
			y_limit2 = m
		yy = y_limit2

	if n1 == 1:
		caption = ax.fill_between(p,c, color="midnightblue", alpha=0.9, label = labelname)
	elif n1 == 2:
		caption = ax.fill_between(p,c, color="crimson", alpha=0.9, label = labelname)

	ax.legend(handles = [caption])
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.set_xticklabels([])
	ax.tick_params(axis='both', bottom=False, which='major', labelsize=8)

	return y_limit


def Generate_annotation_plot(ax, strand, isoforms, exonCountList, exonStartList, exonEndList, startAll, endAll, pos):
	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	rnge = endAll-startAll+1

	ystart = 0
	height = 3
	for i in range(isoforms):
		if i>=15:
			print("15 isoforms of this Gene is plotted.")
			break;
		else:
			ax.hlines(y=(ystart+ystart+height)/2, xmin=startAll, xmax=endAll, linewidth=1, color='skyblue', linestyle = '--')
			
			stList = exonStartList[i]
			enList = exonEndList[i]
			exonCount = int(exonCountList[i])
			for p in range(exonCount):
				ex_s = int(stList[p])
				ex_e = int(enList[p])
				width = int(enList[p]) - int(stList[p]) + 1

				ww = 0.1*width
				if (strand == '+') and (p == exonCount-1):
					width = ex_e - ex_s + 1
					rect1 = patches.Rectangle((ex_s, ystart+1), width, 1, color = 'skyblue', alpha=0.9, fill = True)
					ax.add_patch(rect1)
					rect2 = patches.Rectangle((ex_s, ystart), ww, height, color = 'skyblue', alpha=0.9, fill = True)
					ax.add_patch(rect2)
				elif (strand == '-') and (p==0):
					width = ex_e - ex_s + 1
					rect1 = patches.Rectangle((ex_s, ystart+1), width, 1, color = 'skyblue', alpha=0.9, fill = True)
					ax.add_patch(rect1)
					wstart = ex_s+(0.9*width)
					rect2 = patches.Rectangle((wstart, ystart), ww, height, color = 'skyblue', alpha=0.9, fill = True)
					ax.add_patch(rect2)
				else:
					rect = patches.Rectangle((ex_s,ystart), width, height, color = 'skyblue', fill = True)
					ax.add_patch(rect)

			ystart +=5

	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.set_yticklabels([])
	ax.tick_params(left=False, axis='both', which='major', labelsize=10)

	return

def Plot_Function(pas_flag, region, input1_dir, input2_dir, s1_namelist, s2_namelist, pasSeq1_dir, pasSeq2_dir, p1_namelist, p2_namelist, ann_list, output_dir):
	chrom, geneID, rng = region.split(':')
	start, end = rng.split('-')

	g1_name = input1_dir.split("/")[-1]
	g2_name = input2_dir.split("/")[-1]

	df = pd.DataFrame(ann_list)
	ann_tt = df.loc[df[12]==geneID]
	#if len(ann_tt)==2 or len(ann_tt)==3 or len(ann_tt)==4 or len(ann_tt)==5:
	exonStartList = {}
	exonEndList = {}
	exonCountList = {}

	isoforms = 0
	position = 0
	mini = 500000000
	maxi = 0
	tx_start_list, tx_end_list = [], []
	strand = ""
	
	for a_row in ann_tt.itertuples():
		chrom = a_row[3]
		strand = a_row[4]
		tx_start = int(a_row[5])
		tx_end = int(a_row[6])
		tx_start_list.append(tx_start)
		tx_end_list.append(tx_end)
		cds_start = int(a_row[7])
		cds_end = int(a_row[8])
		exonCount = int(a_row[9])
		exonCountList[isoforms] = exonCount
		exonStartList[isoforms] = ' '.join(a_row[10].split(',')).split()
		exonEndList[isoforms] = ' '.join(a_row[11].split(',')).split()

		isoforms+=1

		if strand == '+':
			if cds_end < mini:
				mini = cds_end
		elif strand == '-':
			if cds_start > maxi:
				maxi = cds_start

	if isoforms > 15:
		print("Cannpt plot genes with more than 15 isoforms")
		sys.exit()

	if strand == '+':
		pos = mini
	else:
		pos = maxi
	all_start = min(np.array(tx_start_list))
	all_end = max(np.array(tx_end_list))

	title = geneID+":"+str(start)+"-"+str(end)

	x_inches = 6.4 
	y_inches = 4.8/5*(len(s1_namelist)*2+1)
	dpi = 100
	fig = plt.figure(1, figsize = (x_inches,y_inches), dpi = dpi, constrained_layout = True)

	number_of_subplots = len(s1_namelist)+len(s2_namelist)+1

	if pas_flag == 0:
		fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1)
		for i in range (0,len(s1_namelist)):
			ax = axes[i]
			if i == 0:
				ax.set_title(title, color = "black", fontsize = 16)
			y_limit = Generate_read_coverate_plot(ax, input1_dir, s1_namelist[i], g1_name, chrom, geneID, all_start, all_end, 1, 0)
		for i in range(len(s1_namelist),number_of_subplots-1):
			j = i - len(s1_namelist)
			y_limit = Generate_read_coverate_plot(axes[i], input2_dir, s2_namelist[j], g2_name, chrom, geneID, all_start, all_end, 2, 0)
	elif pas_flag == 1:
		number_of_subplots = (len(s1_namelist)+len(s2_namelist))*2+1
		if number_of_subplots > 13:
			print("Cannot draw more than 3 samples in each group")
			sys.exit()

		fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1)
		for i in range(0,len(s1_namelist)):
			y_limit = Generate_read_coverate_plot(axes[2*i], input1_dir, s1_namelist[i], g1_name, chrom, geneID, all_start, all_end, 1, 0)
			y_limit2 = Generate_read_coverate_plot(axes[2*i+1], pasSeq1_dir, p1_namelist[i], g1_name, chrom, geneID, all_start, all_end, 1, 1)
			
		for i in range(len(s1_namelist),len(s1_namelist)+len(s2_namelist)):
			print(i)
			j = i - len(s1_namelist)
			y_limit = Generate_read_coverate_plot(axes[2*i], input2_dir, s2_namelist[j], g2_name, chrom, geneID, all_start, all_end, 2, 0)
			y_limit2 = Generate_read_coverate_plot(axes[2*i+1], pasSeq2_dir, p2_namelist[j], g2_name, chrom, geneID, all_start, all_end, 2, 1)
	
	print("Generating annotation plots...")
	ax3 = axes[number_of_subplots-1]
	Generate_annotation_plot(ax3, strand, isoforms, exonCountList, exonStartList, exonEndList, all_start, all_end, pos)
	ax3.set_xlabel('Position', fontsize="14")
	ax3.set_ylabel('Annotation', fontsize="14")

	ax3.spines['top'].set_color('none')
	ax3.spines['bottom'].set_color('none')
	ax3.spines['left'].set_color('none')
	ax3.spines['right'].set_color('none')
	ax3.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

	for i in range(number_of_subplots-1):
		if pas_flag == 0:
			axes[i].set_ylim(0, y_limit*1.1)
		else:
			if i%2==0:
				axes[i].set_ylim(0, y_limit*1.1)
			elif i%2==1:
				axes[i].set_ylim(0, y_limit2*1.1)

	y_limit, y_limit2 = 0, 0
	os.makedirs(output_dir, exist_ok=True)
	plt.savefig(output_dir+title+'.png')
	plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 1000)
	print("Plotted successfully.")


######### Main starts here #################
startTime = time.time()

config = configparser.ConfigParser()
config.read('configuration.ini')

input1_dir = config['INPUT_RNAseq']['input1']
input2_dir = config['INPUT_RNAseq']['input2']
if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]
pasSeq1_dir = config['INPUT_PASseq']['pas1']
pasSeq2_dir = config['INPUT_PASseq']['pas2']
if pasSeq1_dir[:-1] == "/":
	pasSeq1_dir = pasSeq1_dir[:-1]
if pasSeq2_dir[:-1] == "/":
	pasSeq2_dir = pasSeq2_dir[:-1]
output_dir = config['OUTPUT_FOLDER']['output_dir']
if output_dir[-1] != "/":
	output_dir += "/"

os.makedirs(output_dir, exist_ok=True)
inp_annotation = config['ANNOTATION']['annotation']
ref_genome = config['ANNOTATION']['genome']
print("RNA-seq input 1 dir:", input1_dir)
print("RNA-seq input 2 dir:", input2_dir)
print("3'-end-seq input 1 dir:", pasSeq1_dir)
print("3'-end-seq input 2 dir:", pasSeq2_dir)
print("Output Dir:", output_dir) 
print("Annotation:", inp_annotation, ref_genome, "\n\n")

pas_flag = 1
if pasSeq1_dir=='NULL' or pasSeq2_dir == 'NULL':
	pas_flag = 0

s1_namelist = list_dirs(input1_dir)
s2_namelist = list_dirs(input2_dir)
p1_namelist, p2_namelist = "", ""
if pas_flag == 1:
	p1_namelist = list_dirs(pasSeq1_dir)
	p2_namelist = list_dirs(pasSeq2_dir)

print("Loading list of chromosomes from the annotation...")
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

ann_file_reader= open(inp_annotation, "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

region = input("Enter the range: (chr:gene:start-end): ")
Plot_Function(pas_flag, region, input1_dir, input2_dir, s1_namelist, s2_namelist, pasSeq1_dir, pasSeq2_dir, p1_namelist, p2_namelist, ann_list, output_dir)

totalTime = time.time() - startTime
print("Total program time is : ",totalTime)

