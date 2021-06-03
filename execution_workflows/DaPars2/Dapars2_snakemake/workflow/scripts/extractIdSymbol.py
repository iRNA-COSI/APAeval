#!/usr/bin/env python3

import sys

sys.stdout.write("gene_id" + '\t' + "gene_name" + "\n")
pre_attr_id = ""
for line in sys.stdin:
    if not line.startswith("#"):
        attr = dict(item.strip().split(' ') for item in line.split('\t')[8].strip('\n').split(';') if item)
        if attr['gene_id'] != pre_attr_id:
            sys.stdout.write("%s\n" % (attr['gene_id'].strip('\"') + '\t' + attr['gene_name'].strip('\"')))
            pre_attr_id = attr['gene_id']
