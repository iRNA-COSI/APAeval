#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import main as filterPAS
import sys



def main(gt_bed,gtf,min_total_expr_frac, min_frac_site, window_size, out_prefix):
    '''
    '''

    # Read in GTF, subset to exons
    exons = pr.read_gtf(gtf).subset(lambda df: df["Feature"] == "exon")

    # extract terminal exons for each transcript
    exons = exons.assign("exon_number",
                         lambda df: df["exon_number"].astype(float).astype(int))

    t_exons = filterPAS.get_terminal_exons(exons)

    # Merge overlapping terminal exons into union intervals whilst retaining metadata

    # Subset to minimal columns
    t_exons = t_exons[["gene_id", "transcript_id"]]

    print(t_exons)

    # First annotate overlapping intervals (on the same strand) with a common identifier
    # Only terminal exons of the same gene are grouped together in case of overlaps
    # Note: consider requiremnt for common gene_id. Any unintended consequences?
    t_exons = t_exons.cluster(strand=True, by="gene_id")

    # Collapse to union intervals, collapsing contributing transcript IDs to comma-separated string
    m_t_exons = filterPAS.collapse_tes(t_exons,
                                       distinct_cols=["gene_id"],
                                       sep=",")

    m_t_exons = m_t_exons.drop("Cluster")

    print(m_t_exons)




if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],sys.argv[6])
