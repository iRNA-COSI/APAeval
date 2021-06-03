#!/usr/bin/env python

import sys

qapa_results_file=open(sys.argv[1],'r')
qapa_results_file.readline()
qapa_quant_bed=open(sys.argv[2],'w')
header_fields=['chrom','chromStart','chromEnd','name','score','strand']
qapa_quant_bed.write('\t'.join(header_fields)+'\n')
for ln in qapa_results_file:
 ln=ln.strip().split('\t')
 fields=[ln[4],ln[8],ln[9],ln[0],ln[13],ln[7]]
 qapa_quant_bed.write('\t'.join(fields)+'\n')
qapa_quant_bed.close()
