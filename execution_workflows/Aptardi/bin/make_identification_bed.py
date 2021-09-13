#!/usr/bin/env python

import sys

input_gtf_path=sys.argv[1]
input_gtf=open(input_gtf_path,"r")
output_bed=open(input_gtf_path.split(".")[0]+"_identification.bed","w")

for ln in input_gtf:
 if not ln.startswith("#"):
  ln=ln.strip().split("\t")
  if ln[1] == "aptardi": #get the Aptardi identified transcripts only
   attrList=ln[-1].split(";")
   attrDict={}
   for k in attrList:
    p=k.strip().split(" ")
    if len(p) == 2:
     attrDict[p[0]]=p[1].strip('\"')
   name=attrDict["transcript_id"]
   #chrom,chromStart,chromEnd,name,score,strand
   output_bed.write("\t".join([ln[0],ln[3],ln[4],name,".",ln[6]])+"\n")
output_bed.close()
