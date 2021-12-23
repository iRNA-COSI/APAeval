#!/usr/bin/env python3

import sys

sys.stdout.write("#transcript_id" + '\t' + "gene_name" + "\n")
pre_attr_id = ""
for line in sys.stdin:
    if not line.startswith("#"):
        attr = dict(item.strip().split(' ') for item in line.split('\t')[8].strip('\n').split(';') if item)
        if 'transcript_id' in attr and attr['transcript_id'] != pre_attr_id:
            sys.stdout.write("%s\n" % (attr['transcript_id'].strip('\"') + '\t' + attr['gene_name'].strip('\"')))
            pre_attr_id = attr['transcript_id']
