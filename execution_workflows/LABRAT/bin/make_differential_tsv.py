#!/usr/bin/env python

import sys

labrat_results_file=open(sys.argv[1],'r')
results_header=labrat_results_file.readline().strip().split()
labrat_diff_tsv=open(sys.argv[2],'w')
header_fields=['geneID','pval']
labrat_diff_tsv.write('\t'.join(header_fields)+'\n')
geneID_index,pval_index=results_header.index('Gene'),results_header.index('pval')
for ln in labrat_results_file:
 ln=ln.strip().split()
 labrat_diff_tsv.write(ln[geneID_index]+'\t'+ln[pval_index]+'\n')
labrat_diff_tsv
