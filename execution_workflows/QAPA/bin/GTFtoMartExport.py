#!/usr/bin/env python
import sys
gtf_path=sys.argv[1]
gtf=open(gtf_path,'r')
odict={}
for ln in gtf:
 if not ln.startswith('#'):
  ln=ln.strip().split('\t')
  if ln[2] == "transcript":
   attrList=ln[-1].split(";")
   attrDict={}
   for k in attrList:
    p=k.strip().split(" ")
    if len(p) == 2:
     attrDict[p[0]]=p[1].strip('\"')
   tx_id=attrDict["transcript_id"]
   if tx_id not in odict: ##prevent any duplications
    odict[tx_id]=[attrDict["gene_id"].split('.')[0],attrDict["transcript_id"].split('.')[0],
                  attrDict["gene_type"],attrDict["transcript_type"],attrDict["gene_name"]]
outfile=open(gtf_path.split('/')[-1].split('.gtf')[0]+'_martexport.txt','w')
outfile.write('\t'.join(['Gene stable ID','Transcript stable ID','Gene type','Transcript type','Gene name'])+'\n')
for tx_id in odict:
 outfile.write('\t'.join(odict[tx_id])+'\n')
outfile.close()
