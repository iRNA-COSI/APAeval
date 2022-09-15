#!/usr/bin/env python3

import sys

tapa_quant_txt=open(sys.argv[1],'r')
outbed_relative_quant=open(sys.argv[2],'w')
outbed_ident=open(sys.argv[3],'w')

for ln in tapa_quant_txt:
 ln=ln.strip().split('\t')
 g_id,chrom,strand=ln[0].split('NM_')[1],ln[1],ln[2]
 sites,scores=ln[3].split(','),ln[4].split(',')
 scores = [float(score) for score in scores]
 if len(sites) == len(scores):
  for i in range(len(sites)):
   name=g_id+'_'+str(i)
   relative_score=scores[i]/sum(scores)
   if strand == '+':
    chromStart=sites[i]
    chromEnd=str(int(chromStart)+1)
   elif strand == '-':
    chromEnd=sites[i]
    chromStart=str(int(chromEnd)-1)
   outbed_relative_quant.write('\t'.join([chrom,chromStart,chromEnd,name,str(relative_score),strand])+'\n')
   outbed_ident.write('\t'.join([chrom,chromStart,chromEnd,name,'.',strand])+'\n')
outbed_relative_quant.close()
outbed_ident.close()
