#!/usr/bin/env python

import sys

qapa_results_file=open(sys.argv[1],'r')
qapa_results_file.readline()
qapa_quant_bed=open(sys.argv[2],'w')
header_fields=['chrom','chromStart','chromEnd','name','score','strand']
qapa_quant_bed.write('\t'.join(header_fields)+'\n')
for ln in qapa_results_file:
 ln=ln.strip().split('\t')
 strand = ln[7] 
 #if "+" strand gene we want to grab "UTR3.End" as PAS
 if strand == '+': 
  pas_coord = int(ln[9])
 #if "-" strand gene we need "UTR3.Start" as PAS
 else: 
  pas_coord = int(ln[8])
 fields=[ln[4], str(pas_coord-1), str(pas_coord),ln[0],ln[13],ln[7]]
 qapa_quant_bed.write('\t'.join(fields)+'\n')
qapa_quant_bed.close()
