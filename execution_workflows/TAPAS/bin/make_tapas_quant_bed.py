import sys
tapa_quant_txt=open(sys.argv[1],'r')
outbed=open(sys.argv[2],'w')
for ln in tapa_quant_txt:
 ln=ln.strip().split('\t')
 g_id,chrom,strand=ln[0].split('NM_')[1],ln[1],ln[2]
 sites,scores=ln[3].split(','),ln[4].split(',')
 if len(sites) == len(scores):
  for i in range(len(sites)):
   name=g_id+'_'+str(i)
   score=scores[i]
   if strand == '+':
    chromStart=sites[i]
    chromEnd=str(int(chromStart)+1)
   elif strand == '-':
    chromEnd=sites[i]
    chromStart=str(int(chromEnd)+1)
   outbed.write('\t'.join([chrom,chromStart,chromEnd,name,score,strand])+'\n')
outbed.close()
