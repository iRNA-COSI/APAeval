#!/usr/bin/env python3

import sys
bed_path,gtf_path=sys.argv[1],sys.argv[2]
def readGTF(gtf_path):
 transcript2gene={}
 gtf=open(gtf_path,"r")
 for ln in gtf:
  if not ln.startswith("#"):
   attrList=ln.split("\t")[-1].split(";")
   attrDict={}
   for k in attrList:
    p=k.strip().split(" ")
    if len(p) == 2:
     attrDict[p[0]]=p[1].strip('\"')
   if "transcript_id" in attrDict:
    transcript2gene[attrDict["transcript_id"]]=attrDict["gene_id"]
 return transcript2gene
transcript2gene=readGTF(gtf_path)

tapas_ref_path="tapas_ref.txt"
bed=open(bed_path,"r")
outfile=open(tapas_ref_path,"w")
for ln in bed:
 ln=ln.split("\t")
 transcript="NM_"+ln[3]
 gene="NM_"+transcript2gene[ln[3]]
 organized=[gene,transcript,ln[0],ln[5],ln[1],ln[2],ln[6],ln[7],ln[9],ln[10],ln[11]]
 outfile.write("\t".join(organized))
outfile.close()
