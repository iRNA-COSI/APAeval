#!/usr/bin/env python

import sys

qapa_results_file=open(sys.argv[1],'r')
qapa_results_file.readline()
qapa_quant_bed_ppau=open(sys.argv[2],'w')
qapa_quant_bed_tpm=open(sys.argv[3],'w')

for ln in qapa_results_file:
 ln=ln.strip().split('\t')
 # if the ppau value is NA, don't calculate ppau fraction
 if ln[12] == "NA":
  ppau_fraction = ln[12]
 # if ppau is not NA, calculate ppau fraction by dividing by 100
 else:
  ppau_fraction=float(ln[12])/100.0
 if ln[7] == '+':
  chromStart,chromEnd=str(int(ln[9])-1),ln[9]
 elif ln[7] == '-':
  chromStart,chromEnd=ln[8],str(int(ln[8])+1)
 ##write to the ppau bed file
 ppau_fields=[ln[4],chromStart,chromEnd,ln[0],str(ppau_fraction),ln[7]]
 qapa_quant_bed_ppau.write('\t'.join(ppau_fields)+'\n')
 ##write to the tmp bed file
 tpm_fields=[ln[4],chromStart,chromEnd,ln[0],ln[13],ln[7]]
 qapa_quant_bed_tpm.write('\t'.join(tpm_fields)+'\n')
qapa_quant_bed_ppau.close()
qapa_quant_bed_tpm.close()
