#!/usr/bin/env python

import sys

qapa_results_file=open(sys.argv[1],'r')
qapa_results_file.readline()
qapa_quant_bed=open(sys.argv[2],'w')

for ln in qapa_results_file:
 ln=ln.strip().split('\t')
 ppau_fraction=float(ln[12])/100.0
 if ln[7] == '+':
  chromStart,chromEnd=str(int(ln[9])-1),ln[9]
 elif ln[7] == '-':
  chromStart,chromEnd=ln[8],str(int(ln[8])+1)
 fields=[ln[4],chromStart,chromEnd,ln[0],str(ppau_fraction),ln[7]]
 qapa_quant_bed.write('\t'.join(fields)+'\n')
qapa_quant_bed.close()
