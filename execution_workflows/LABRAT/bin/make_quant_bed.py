#!/usr/bin/env python

import sys

labrat_quant_file=open(sys.argv[1],'r')
labrat_quant_file.readline()
labrat_quant_bed=open(sys.argv[2],'w')
header_fields=['chrom','chromStart','chromEnd','name','score','strand']
labrat_quant_bed.write('\t'.join(header_fields)+'\n')
gff=open(sys.argv[3],'r')
def getChrStrand(gff):
 dict={}
 for ln in gff:
  if not ln.startswith('#'): 
   ln=ln.strip().split()
   if ln[2] == 'transcript':
    chrom,strand,tx_id=ln[0],ln[6],ln[-1].split('transcript_id=')[1].split('.')[0] ##salmon quant does not include version
    if tx_id not in dict:
     dict[tx_id]={'chrom':chrom,'strand':strand}
 return dict

annotation_dict=getChrStrand(gff)
for ln in labrat_quant_file:
 ln=ln.strip().split('\t')
 tx_id,TPM=ln[0],ln[3]
 fields=[annotation_dict[tx_id]['chrom'],'NA','NA',tx_id,TPM,annotation_dict[tx_id]['strand']]
 labrat_quant_bed.write('\t'.join(fields)+'\n')
labrat_quant_bed.close()
