#!/usr/bin/bash

#########################################################
### Generate results for identification & quantification
#########################################################

SAMPLE=${1} # sample id
ASSEMDIR=${2} # the assembly dir in isoscm result, isoscm/ by default


# identification
cat ${ASSEMDIR}/tmp/${SAMPLE}.cp.filtered.gtf | awk -vOFS="\t" \
          '{print $1,$4,$5,"site"NR,$8,$7}' > ${SAMPLE}_IsoSCM_01.bed

# quantification
cat ${ASSEMDIR}/${SAMPLE}.coverage.gtf | awk '$9 ~ /coverage/' | \
    sed -E 's/(.*coverage ")([^"]*)(.*)/\1\2\3\t\2/g' > tmp.gtf

awk -F"\t" -vOFS="\t" 'NR==FNR{bed[$1"_"$2]=$0;next} { for (b in bed) \
    {split(b,pos,"_"); if (pos[1]==$1 && pos[2]>=$4 && pos[2]<=$5) \
    {print bed[b]"\t"$10; delete bed[b]}} }' ${SAMPLE}_IsoSCM_01.bed tmp.gtf | \
    awk -vOFS="\t" '{print $1,$2,$3,$4,$7,$6}' > ${SAMPLE}_IsoSCM_02.bed

