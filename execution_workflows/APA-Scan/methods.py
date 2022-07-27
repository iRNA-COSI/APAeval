import csv
import time
import sys
import math
from operator import itemgetter
import pandas as pd
from Bio import SeqIO
import re
import bisect
from bisect import bisect_left
from scipy.stats import chisquare
import peakutils
import numpy as np
import os
import glob, os #####

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()

def bi_contains(lst, item):
    return bisect_left(lst, item)

def getFinalTargetRegion(inputList):
	n = len(inputList)
	inputList.sort(key = itemgetter(1), reverse = True)
	st = Stack()
	st.push(inputList[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputList[i][0] < stacktop[0] and inputList[i][1] < stacktop[0]:
			st.push(inputList[i])
		elif inputList[i][0] < stacktop[0] and inputList[i][1] > stacktop[0] and inputList[i][1] < stacktop[1]:
			st.pop()
			st.push((inputList[i][0], stacktop[1]))
		elif inputList[i][1] == stacktop[1] and inputList[i][0] < stacktop[0]:
			st.pop()
			st.push(inputList[i])

	newList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		newList.append(stacktop)
		st.pop()

	return newList


def findAllOccurance(p, s):
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def CountReadCoverage(chrom, start, end, bam_list, position_row):
	totalCount = 0

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if pos2 >= len(bam_list) or int(bam_list[pos2][1]) != end:
		pos2 = pos2 - 1

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		for t in range(pos1, pos2+1):
			read = int(bam_list[t][2])
			totalCount += read

	return totalCount

def read_bamfiles(input_dir, sample, chrom):
	ss = time.time()
	bam_file_reader = open(input_dir+'/'+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]
	return bam_list, position_row

def Generate_coverage(chrom, start, end, pos, strand, bam_list, position_row):
	if strand == '+':
		length = pos-start
		targetLength = end-pos
		RC = CountReadCoverage(chrom, start, pos, bam_list, position_row)
		targetRC = CountReadCoverage(chrom, pos+1, end, bam_list, position_row)
	else:
		targetLength = pos-start
		length = end-pos
		RC = CountReadCoverage(chrom, pos+1, end, bam_list, position_row)
		targetRC = CountReadCoverage(chrom, start, pos, bam_list, position_row)

	n = targetRC/targetLength
	N = RC/length
	return n, N

def with_PAS_signal(chromosomes, input1_dir, input2_dir, s1_namelist, s2_namelist, g1, g2, output_dir, result_filename):
	ss = time.time()
	len1 = len(s1_namelist)
	len2 = len(s2_namelist)

	read1 = open(output_dir+"Signal_positions.csv", "r")
	reader1 = csv.reader(read1, delimiter="\t")
	pasList = list(reader1)
	df_p = pd.DataFrame(pasList)

	with open(result_filename+".csv",'w') as f: ### HAD TO REMOVE INITIAL "output_dir+" FROM STRING CUS OUTPUT STRING IS INCLUDED IN "result_filename" ~ EM
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'n1: '+g1, 'n2: '+g2, 'N1: '+g1, 'N2: '+g2])

		position_row = []
		for chrom in chromosomes:
			s1_bam_list = {}
			s2_bam_list = {}
			s1_position_row = {}
			s2_position_row = {}
			for sample1 in s1_namelist:
				bam_list, position_row = read_bamfiles(input1_dir, sample1, chrom)
				s1_bam_list[sample1] = bam_list
				s1_position_row[sample1] = position_row
			for sample2 in s2_namelist:
				bam_list, position_row = read_bamfiles(input2_dir, sample2, chrom)
				s2_bam_list[sample2] = bam_list
				s2_position_row[sample2] = position_row

			selected_rows = df_p.loc[df_p[0]==chrom]
			for row in selected_rows.itertuples():
				geneID = row[2].strip()
				strand = row[3]
				start = int(row[4])
				end = int(row[5])
				pasPositionsList = (str(row[6]).strip('[ ]')).split(',')

				signi_p = 1.1
				ratio_diff = 'Nan'
				signi_ratio_diff = 'Nan'
				abs_ratio_diff = 'Nan'
				n1_f = 0
				N1_f = 0
				n2_f = 0
				N2_f = 0
				signi_pos = 0
				flag = 0

				for pos in pasPositionsList:
					pos = int(pos.strip())
					length = end - start + 1
					targetLength = length

					s1_n, s1_N, s2_n, s2_N = 0, 0, 0, 0
					for sample1 in s1_namelist:
						n, N = Generate_coverage(chrom, start, end, pos, strand, s1_bam_list[sample1], s1_position_row[sample1])
						s1_n += n
						s1_N += N
					for sample2 in s2_namelist:
						n, N = Generate_coverage(chrom, start, end, pos, strand, s2_bam_list[sample2], s2_position_row[sample2])
						s2_n += n
						s2_N += N

					n1 = s1_n/len1
					N1 = s1_N/len1
					n2 = s2_n/len2
					N2 = s2_N/len2

					ratio_diff = 0

					N1 = N1 + n1
					N2 = N2 + n2

					if N1>0 and N2>0:
						ratio_diff = (n1/N1) - (n2/N2)
						P0 = (n1+n2)/(N1+N2)
						n10 = N1 * P0
						n20 = N2 * P0
						exp = [n10, N1-n10, n20, N2-n20]
						if 0 not in exp:
							flag = 1
							res, p_value = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
							if p_value < signi_p:
								signi_p = p_value
								signi_ratio_diff = ratio_diff
								abs_ratio_diff = abs(signi_ratio_diff)
								n1_f = n1
								N1_f = N1
								n2_f = n2
								N2_f = N2
								signi_pos = pos

				if flag == 1:
					writer.writerow([chrom, geneID, strand, start, end, signi_pos, signi_p, signi_ratio_diff, abs_ratio_diff, n1_f, n2_f, N1_f-n1_f, N2_f-n2_f])
			print(chrom, "completed in ", chrom, time.time() - ss, "sec")
		f.close()
	os.remove(output_dir+"Signal_positions.csv") ### HAD TO DELTE "r" OPTION HERE AS :- "TypeError: Function takes at most 1 positional arguments (2 given)" ~ "EM"
	return

def makeChromDict(chromosomes, inp_annotation):
	chromDict = {}
	geneFlag = []
	with open(inp_annotation, 'r') as f:
	    reader = csv.reader(f, dialect='excel', delimiter='\t')
	    headers = next(f)
	    readerList = list(reader)
	    for chrom in chromosomes:
	    	geneList = []
	    	df = pd.DataFrame(readerList)
	    	rowsOfChr = df.loc[df[2] == chrom]
	    	for row in rowsOfChr.itertuples():
	    		geneId = row[13].strip()
	    		strand = row[4]
	    		if geneId not in geneFlag:
	    			geneList.append((geneId, strand))
	    			geneFlag.append(geneId)

	    	chromDict[chrom] = geneList

	return chromDict


def getFinalTargetRegion(inputList):
	n = len(inputList)
	inputList.sort(key = itemgetter(1), reverse = True)
	st = Stack()
	st.push(inputList[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputList[i][0] < stacktop[0] and inputList[i][1] < stacktop[0]:
			st.push(inputList[i])
		elif inputList[i][0] < stacktop[0] and inputList[i][1] > stacktop[0] and inputList[i][1] < stacktop[1]:
			st.pop()
			st.push((inputList[i][0], stacktop[1]))
		elif inputList[i][1] == stacktop[1] and inputList[i][0] < stacktop[0]:
			st.pop()
			st.push(inputList[i])

	newList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		newList.append(stacktop)
		st.pop()

	return newList


def makeSplittedList(position_row, bam_list, start, end):
	p = []
	r = []

	pos1 = bi_contains(position_row, start)
	pos2 = bi_contains(position_row, end)

	if pos2 >= len(bam_list) or int(bam_list[pos2][1]) != end:
		pos2 = pos2 - 1

	if(pos1 < len(bam_list) and pos2 < len(bam_list)):
		for t in range(pos1, pos2+1):
			p.append(int(bam_list[t][1]))
			r.append(int(bam_list[t][2]))

	return (p,r)

def calculatePeaksWithVallys(peakAreaPos, peakAreaRead, indexes, allStartPoint, allEndPoint):
	peakRangeList = []
	peakStartPoint = allStartPoint
	prevPeakValue = 0
	indexLen = len(indexes)

	# because we will compare two peaks together to consider the valley in between, we wont take the last peak index in the for loop
	for k in range(indexLen-1):
		selectedRange = peakAreaRead[indexes[k]:indexes[k+1]+1]
		valley = min(selectedRange)
		valleyPos = selectedRange.index(valley)
		actualValleyPos = peakAreaPos[indexes[k]+valleyPos]

		if valley < 0.3*min(peakAreaRead[indexes[k]],peakAreaRead[indexes[k+1]]):
			peakEndPoint = actualValleyPos
			peakRangeList.append((peakAreaPos[peakAreaRead.index(max(prevPeakValue,peakAreaRead[indexes[k]]))], peakStartPoint, peakEndPoint))
			prevPeakValue = 0
			peakStartPoint = actualValleyPos + 1
		else:
			prevPeakValue = max(prevPeakValue, peakAreaRead[indexes[k]])


	ind = [i for i, value in enumerate(peakAreaRead) if value == max(prevPeakValue, peakAreaRead[indexes[indexLen-1]]) and peakAreaPos[i] >= peakStartPoint]

	# because ind is a list with indexes of the maxpeak value, for that we will take ind[0] (ind[0] or ind[1] are the indices of the same value)
	peakRangeList.append((peakAreaPos[ind[0]], peakStartPoint, allEndPoint))
	return peakRangeList


def findPeakPosition(p, r):
	listOfPeakPositions = []
	listOfPeakRange = []
	flag = 0
	for i in range(len(p)):
		if int(r[i])>1:
			if flag == 0:
				peakStartPoint = p[i]
				peakAreaPos = []
				peakAreaRead = []
				flag = 1
			if flag == 1:
				peakAreaPos.append(p[i])
				peakAreaRead.append(r[i])

		if (int(r[i])==1 and flag == 1) or (i == len(p)-1 and flag == 1):
			if i == len(p)-1:
				peakEndPoint = p[i]
			else:
				peakEndPoint = p[i-1]

			np_peakAreaPos = np.array(peakAreaPos)
			np_peakAreaRead = np.array(peakAreaRead)

			indexes = peakutils.indexes(np_peakAreaRead, thres=6, min_dist=35, thres_abs = True)

			numberOfIndex = len(indexes)

			if numberOfIndex == 1:
				listOfPeakRange.append((peakAreaPos[indexes[0]], peakStartPoint, peakEndPoint))

			elif numberOfIndex > 1:
				peakRangeList = calculatePeaksWithVallys(peakAreaPos, peakAreaRead, indexes, peakStartPoint, peakEndPoint)

				listOfPeakRange.extend(peakRangeList)

			flag = 0

	return listOfPeakRange

def mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand):
	len1 = len(listOfPeakRange1)
	len2 = len(listOfPeakRange2)

	cleavageSites = []
	peakFlag = []
	minVal = 10000

	if len1 == 0:
		for (peak, st, en) in listOfPeakRange2:
			if strand == '+':
				cleavageSites.append(en)
			else:
				cleavageSites.append(st)
		return cleavageSites
	elif len2 == 0:
		for (peak, st, en) in listOfPeakRange1:
			if strand == '+':
				cleavageSites.append(en)
			else:
				cleavageSites.append(st)
		return cleavageSites

	if len1<len2:
		row = len1
		col = len2
		firstList = listOfPeakRange1
		secondList = listOfPeakRange2
	else:
		row = len2
		col = len1
		firstList = listOfPeakRange2
		secondList = listOfPeakRange1

	peakFlag1 = {}
	peakFlag2 = {}
	distanceToSortList = []

	for i in range(row):
		for j in range(col):
			if strand == '+':
				distance = abs(int(firstList[i][2])-int(secondList[j][2]))
				distanceToSortList.append((int(firstList[i][2]), int(secondList[j][2]), distance))
				peakFlag1[firstList[i][2]] = 0
				peakFlag2[secondList[j][2]] = 0
			else:
				distance = abs(int(firstList[i][1])-int(secondList[j][1]))
				distanceToSortList.append((int(firstList[i][1]), int(secondList[j][1]), distance))
				peakFlag1[firstList[i][1]] = 0
				peakFlag2[secondList[j][1]] = 0

	distanceToSortList.sort(key = itemgetter(2), reverse = False)

	for (pos1, pos2, dist) in distanceToSortList:
		if peakFlag1[pos1] == 0 and peakFlag2[pos2] == 0:
			averageOfTwoPeaks = math.ceil((pos1+pos2)/2)
			cleavageSites.append(averageOfTwoPeaks)
			peakFlag1[pos1] = 1
			peakFlag2[pos2] = 1
		elif peakFlag1[pos1] != 0 and peakFlag2[pos2] == 0:
			cleavageSites.append(pos2)
			peakFlag2[pos2] = 1
		elif peakFlag1[pos1] == 0 and peakFlag2[pos2] != 0:
			cleavageSites.append(pos1)
			peakFlag1[pos1] = 1

	return cleavageSites


def Generate_withPasSeqData(filename, chromosomes, inp_annotation, p1_dir, p2_dir, p1_name, p2_name, output_dir):
	chromDict = makeChromDict(chromosomes, inp_annotation)

	with open(inp_annotation, 'r') as f:
	    reader = csv.reader(f, dialect='excel', delimiter='\t')
	    headers = next(f)
	    readerList = list(reader)
	    df = pd.DataFrame(readerList)

	    with open(filename,'w') as g:
	    	writer = csv.writer(g, delimiter='\t')
	    	writer.writerow(['Chrom', 'Gene Name', 'Strand', 'Start', 'End',  'Cleavage Sites'])

	    	for chrom in chromosomes:
	    		tt = time.time()
	    		bam_file_reader1 = open(p1_dir+"/"+p1_name+"/"+chrom+".txt", "rt")
	    		bam_read1 = csv.reader(bam_file_reader1, delimiter="\t")
	    		bam_list1 = list(bam_read1)
	    		position_row1 = [int(bam_list1[i][1]) for i in range(len(bam_list1))]

	    		bam_file_reader2 = open(p2_dir+"/"+p2_name+"/"+chrom+".txt", "rt")
	    		bam_read2 = csv.reader(bam_file_reader2, delimiter="\t")
	    		bam_list2 = list(bam_read2)
	    		position_row2 = [int(bam_list2[i][1]) for i in range(len(bam_list2))]

	    		geneList = chromDict[chrom]
	    		for (geneId, strand) in geneList:
	    			rowsOfChr = df.loc[(df[2] == chrom) & (df[12] == geneId)]
	    			targetList = []
	    			for row in rowsOfChr.itertuples():
	    				if str(row[2].strip()).startswith('NM_') or str(row[2].strip()).startswith('NR_'):
		    				exonCount = int(row[9])
		    				exonStartList = row[10].split(',')
		    				exonEndList = row[11].split(',')
		    				txStart = int(row[5])
		    				txEnd = int(row[6])
		    				cdsStart = int(row[7])
		    				cdsEnd = int(row[8])

		    				if cdsStart < cdsEnd:
			    				if strand == '+':
			    					for i in range(exonCount):
			    						exonStart = int(exonStartList[i])
			    						exonEnd = int(exonEndList[i])
			    						if(exonStart<=cdsEnd and cdsEnd<=exonEnd):
			    							if (exonStart, txEnd) not in targetList:
			    								targetList.append((exonStart, txEnd))

			    				elif strand == '-':
			    					for i in range(exonCount):
			    						exonStart = int(exonStartList[i])
			    						exonEnd = int(exonEndList[i])
			    						if(exonStart<=cdsStart and cdsStart<=exonEnd):
			    							if (txStart, exonEnd) not in targetList:
			    								targetList.append((txStart, exonEnd))

		    		if len(targetList) > 0:
		    			revisedTargetList = getFinalTargetRegion(targetList)

		    			for (st, en) in revisedTargetList:
		    				(p1, r1) = makeSplittedList(position_row1, bam_list1, st, en)
		    				(p2, r2) = makeSplittedList(position_row2, bam_list2, st, en)

		    				listOfPeakRange1 = findPeakPosition(p1, r1)
		    				listOfPeakRange2 = findPeakPosition(p2, r2)

		    				if len(listOfPeakRange1) > 1 or len(listOfPeakRange2) > 1:

		    					cleavageSites = mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand)
		    					writer.writerow([chrom, geneId, strand, st, en, cleavageSites])

		    	print("Completed processing",chrom, ". Total time:", round(time.time() - tt, 2), "sec.")

	    g.close()
	f.close()


def Quantification(chromosomes, s1_dir, s2_dir, g1_name, g2_name, filename, result_filename):
	ss = time.time()

	with open(filename, 'r') as h:
		reader = csv.reader(h, dialect='excel', delimiter='\t')
		headers = next(h)
		pasReadList = list(reader)

		with open(result_filename+".csv",'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['Chrom', 'Gene Name', 'strand', 'Start', 'End', 'Position', 'p-value', 'Ratio Difference', 'Absolute ratio difference', 'Length (before pos)', 'Read count (before pos)'+g1_name, 'Read count (before_pos)'+g2_name, 'Length (After_pos)', 'Read Count (After pos)'+g1_name, 'Read Count (After pos)'+g2_name, g1_name+':n1', g1_name+':N1', g2_name+':n2', g2_name+':N2'])

			position_row = []
			for chrom in chromosomes:
				bam_file_reader1 = open(s1_dir+"/"+g1_name+'/'+chrom+".txt", "rt")
				bam_read1 = csv.reader(bam_file_reader1, delimiter="\t")
				bam_list1 = list(bam_read1)
				position_row1 = [int(bam_list1[i][1]) for i in range(len(bam_list1))]

				bam_file_reader2 = open(s2_dir+"/"+g2_name+'/'+chrom+".txt", "rt")
				bam_read2 = csv.reader(bam_file_reader2, delimiter="\t")
				bam_list2 = list(bam_read2)
				position_row2 = [int(bam_list2[i][1]) for i in range(len(bam_list2))]

				df = pd.DataFrame(pasReadList)
				target_tt = df.loc[df[0]==chrom]
				for t_row in target_tt.itertuples():
					geneID = t_row[2].strip()
					strand = t_row[3]
					start = int(t_row[4])
					end = int(t_row[5])
					pasPositionsList = (str(t_row[6]).strip('[ ]')).split(',')

					signi_p = 1.1
					ratio_diff = 'Nan'
					signi_ratio_diff = 'Nan'
					abs_ratio_diff = 'Nan'
					n1_f = 0
					N1_f = 0
					n2_f = 0
					N2_f = 0
					targetRC1_f = 0
					targetRC2_f = 0
					RC1_f = 0
					RC2_f = 0
					length_f = 0
					targetLength_f = 0
					signi_pos = 0
					flag = 0
					for pos in pasPositionsList:
						pos = int(pos.strip())
						length = end - start + 1
						targetLength = length

						if (strand == '+' and (pos-start)<=(length*0.85)) or (strand == '-' and (pos-start)>=(length*0.15)):
							if strand == '+':
								length = pos-start
								targetLength = end - pos
								RC1 = CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
								RC2 = CountReadCoverage(chrom, start, pos, bam_list2, position_row2)
								targetRC1 = CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
								targetRC2 = CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
							else:
								targetLength = pos-start
								length = end - pos
								RC1 = CountReadCoverage(chrom, pos+1, end, bam_list1, position_row1)
								RC2 = CountReadCoverage(chrom, pos+1, end, bam_list2, position_row2)
								targetRC1 = CountReadCoverage(chrom, start, pos, bam_list1, position_row1)
								targetRC2 = CountReadCoverage(chrom, start, pos, bam_list2, position_row2)

							n1 = targetRC1/targetLength
							n2 = targetRC2/targetLength
							N1 = RC1/length
							N2 = RC2/length

							if N1!=0 and N2!=0:
								ratio_diff = (n1/N1) - (n2/N2)
							else:
								ratio_diff = 0

							N1 = N1 + n1
							N2 = N2 + n2

							if N1!= 0 and N2!=0:
								P0 = (n1+n2)/(N1+N2)
								n10 = N1 * P0
								n20 = N2 * P0
								exp = [n10, N1-n10, n20, N2-n20]
								if 0 not in exp:
									flag = 1
									res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
									if res[1] < signi_p:
										signi_p = res[1]
										signi_ratio_diff = ratio_diff
										abs_ratio_diff = abs(signi_ratio_diff)
										n1_f = n1
										N1_f = N1
										n2_f = n2
										N2_f = N2
										RC1_f = RC1
										RC2_f = RC2
										targetRC1_f = targetRC1
										targetRC2_f = targetRC2
										signi_pos = pos
										length_f = length
										targetLength_f = targetLength

					if flag == 1:
						writer.writerow([chrom, geneID, strand, start, end, signi_pos, signi_p, signi_ratio_diff, abs_ratio_diff, length_f, RC1_f, RC2_f, targetLength_f, targetRC1_f, targetRC2_f, n1_f, N1_f-n1_f, n2_f, N2_f-n2_f])
			f.close()
		h.close()
	print("End of quantification")

######################################################################################
def with_pasSeq(chromosomes, chromDict, pasSeq_dir1, pasSeq_dir2, p1_namelist, p2_namelist, g1_name, g2_name, inp_annotation, output_dir):
	print()


def Get_Peak_Positions(chromosomes, chromDict, pasSeq_dir, pas_sample, inp_annotation, output_dir):
	with open(inp_annotation, 'r') as f:
	    reader = csv.reader(f, dialect='excel', delimiter='\t')
	    headers = next(f)
	    readerList = list(reader)
	    df = pd.DataFrame(readerList)

	    filename = output_dir+pas_sample+'_peak_positions.csv'
	    with open(filename,'w') as g:
	    	writer = csv.writer(g, delimiter='\t')
	    	writer.writerow(['Chrom', 'Gene Name', 'Strand', 'Start', 'End',  'Cleavage Sites'])

	    	for chrom in chromosomes:
	    		tt = time.time()
	    		bam_file_reader1 = open(pasSeq_dir+"/"+pas_sample+'/'+chrom+".txt", "rt")
	    		bam_read1 = csv.reader(bam_file_reader1, delimiter="\t")
	    		bam_list1 = list(bam_read1)
	    		position_row1 = [int(bam_list1[i][1]) for i in range(len(bam_list1))]

	    		geneList = chromDict[chrom]
	    		for (geneId, strand) in geneList:
	    			rowsOfChr = df.loc[(df[2] == chrom) & (df[12] == geneId)]
	    			targetList = []
	    			for row in rowsOfChr.itertuples():
	    				if str(row[2].strip()).startswith('NM_') or str(row[2].strip()).startswith('NR_'):
		    				exonCount = int(row[9])
		    				exonStartList = row[10].split(',')
		    				exonEndList = row[11].split(',')
		    				txStart = int(row[5])
		    				txEnd = int(row[6])
		    				cdsStart = int(row[7])
		    				cdsEnd = int(row[8])

		    				if cdsStart < cdsEnd:
			    				if strand == '+':
			    					for i in range(exonCount):
			    						exonStart = int(exonStartList[i])
			    						exonEnd = int(exonEndList[i])
			    						if(exonStart<=cdsEnd and cdsEnd<=exonEnd):
			    							if (exonStart, txEnd) not in targetList:
			    								targetList.append((exonStart, txEnd))

			    				elif strand == '-':
			    					for i in range(exonCount):
			    						exonStart = int(exonStartList[i])
			    						exonEnd = int(exonEndList[i])
			    						if(exonStart<=cdsStart and cdsStart<=exonEnd):
			    							if (txStart, exonEnd) not in targetList:
			    								targetList.append((txStart, exonEnd))

		    		if len(targetList) > 0:
		    			revisedTargetList = getFinalTargetRegion(targetList)

		    			for (st, en) in revisedTargetList:
		    				(p1, r1) = makeSplittedList(position_row1, bam_list1, st, en)
		    				(p2, r2) = makeSplittedList(position_row2, bam_list2, st, en)

		    				listOfPeakRange1 = findPeakPosition(p1, r1)
		    				listOfPeakRange2 = findPeakPosition(p2, r2)

		    				if len(listOfPeakRange1) > 1 or len(listOfPeakRange2) > 1:

		    					cleavageSites = mergePeaksFromBothSamples(listOfPeakRange1, listOfPeakRange2, strand)
		    					writer.writerow([chrom, geneId, strand, st, en, cleavageSites])

		    	print("Completed processing",chrom, ". Total time:", round(time.time() - tt, 2), "sec.")

	    g.close()
	f.close()


def Get_Signal_Positions(chromosomes, chromDict, inp_annotation, ref_genome, output_dir):
	with open(inp_annotation, 'r') as f:
	    reader = csv.reader(f, dialect='excel', delimiter='\t')
	    headers = next(f)
	    readerList = list(reader)
	    df = pd.DataFrame(readerList)

	    filename = output_dir+'Signal_positions.csv'
	    with open(filename,'w') as g:
	    	writer = csv.writer(g, delimiter='\t')
	    	writer.writerow(['Chrom', 'Gene Name', 'Strand', 'Start', 'End', 'PAS positions'])

	    	fasta_sequences = SeqIO.parse(open(ref_genome),'fasta')
	    	for fasta in fasta_sequences:
	    		chrom, sequence = fasta.id, str(fasta.seq)
	    		if chrom in chromosomes:
	    			geneList = chromDict[chrom]
	    			tt = time.time()
	    			for (geneId, strand) in geneList:
	    				rowsOfChr = df.loc[(df[2] == chrom) & (df[12] == geneId)]
	    				targetList = []
	    				for row in rowsOfChr.itertuples():
	    					if str(row[2].strip()).startswith('NM_') or str(row[2].strip()).startswith('NR_'):
		    					exonCount = int(row[9])
		    					exonStartList = row[10].split(',')
		    					exonEndList = row[11].split(',')
		    					exonStartList_int = [int(exonStartList[i]) for i in range(exonCount)]
		    					exonEndList_int = [int(exonEndList[i]) for i in range(exonCount)]
		    					txStart = int(row[5])
		    					txEnd = int(row[6])
		    					cdsStart = int(row[7])
		    					cdsEnd = int(row[8])
		    					if cdsStart < cdsEnd:
			    					if strand == '+':
			    						regStart = cdsEnd
			    						regEnd = txEnd
			    						for i in range(exonCount):
			    							exonStart = int(exonStartList[i])
			    							exonEnd = int(exonEndList[i])

			    							if(exonStart<=regStart and regStart<=exonEnd):
			    								if (exonStart, regEnd) not in targetList:
			    									targetList.append((exonStart, regEnd))
			    									break

			    					elif strand == '-':
			    						regStart = txStart
			    						regEnd = cdsStart
			    						for i in range(exonCount):
			    							exonStart = int(exonStartList[i])
			    							exonEnd = int(exonEndList[i])
			    							if(exonStart<=regEnd and regEnd<=exonEnd):
			    								if (regStart, exonEnd) not in targetList:
			    									targetList.append((regStart, exonEnd))
			    									break

		    			if len(targetList) > 0:
		    				revisedTargetList = getFinalTargetRegion(targetList)

		    				if strand == '+':
		    					for (st, en) in revisedTargetList:
		    						length = en - st + 1
		    						seq = sequence[st-2: en-2].upper()
		    						listPlus = []
		    						for i in findAllOccurance('AATAAA', seq):
		    							pos = i+5
		    							if pos<=(length*0.85):
			    							listPlus.append(st+pos)
		    						for i in findAllOccurance('ATTAAA', seq):
		    							pos = i+5
		    							if pos<=(length*0.85):
			    							listPlus.append(st+pos)

		    						if len(listPlus)>0:
		    							writer.writerow([chrom, geneId, strand, st, en, list(set(listPlus))])
		    				elif strand == '-':
		    					for (st, en) in revisedTargetList:
		    						length = en - st + 1
		    						seq = sequence[st-2: en-2].upper()
		    						listMinus = []
		    						for pos in findAllOccurance('TTTATT', seq):
		    							if pos>=(length*0.15):
		    								listMinus.append(st+pos)
		    						for pos in findAllOccurance('TTTAAT', seq):
		    							if pos>=(length*0.15):
			    							listMinus.append(st+pos)

			    					if len(listMinus)>0:
			    						writer.writerow([chrom, geneId, strand, st, en, list(set(listMinus))])

	    			print("Completed processing",chrom, ". Total time:", round(time.time() - tt, 2), "sec.")

	    g.close()
	f.close()
