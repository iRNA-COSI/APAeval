#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import main as filterPAS
import sys


def group_sum(df, group_col, value_col, out_col):
    '''
    Calculate a group-wise sum
    df: pd.DataFrame
    group_col: column containing group identifier
    value_col: column label of value to sum group-wise
    returns: 2 column df of group_col | out_col where out_col contains the sum of value_col values for every group
    '''

    sum_df = (df.groupby(group_col)
              [value_col]
              .sum()
              .reset_index()
              .rename(columns={value_col: out_col})
              )

    return sum_df


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

    # print(t_exons)

    # First annotate overlapping intervals (on the same strand) with a common identifier
    # Only terminal exons of the same gene are grouped together in case of overlaps
    # Note: consider requiremnt for common gene_id. Any unintended consequences?
    t_exons = t_exons.cluster(strand=True, by="gene_id")

    # Collapse to union intervals, collapsing contributing transcript IDs to comma-separated string
    m_t_exons = filterPAS.collapse_tes(t_exons,
                                       distinct_cols=["gene_id"],
                                       sep=",")

    m_t_exons = m_t_exons.drop("Cluster")

    # print(m_t_exons)

    # Construct the minimal terminal exon ID - gene_id|<transcript_id1,transcript_id2>

    m_t_exons = filterPAS.assign_id(m_t_exons,
                                    cols_to_cat=["gene_id", "transcript_id"],
                                    sep_char="|",
                                    out_col="te_id_min")

    m_t_exons = m_t_exons.drop(["gene_id", "transcript_id"])

    # print(m_t_exons)

    # read in ground truth polyA sites
    pas = pr.read_bed(gt_bed)

    # assign an id column (pas_id) so can track number of PAS
    # Rename score column so clear it corresponds to TPM
    pas = filterPAS.assign_id(pas)
    pas = pas.apply(lambda df: df.rename(columns={"Score": "tpm"}))

    # Find TEs with at least two overlapping PAS
    print(f"Number of union terminal exons - {m_t_exons.te_id_min.nunique()}")
    m_t_exons = m_t_exons.count_overlaps(pas, strandedness="same")
    m_t_exons = m_t_exons.subset(lambda df: df["NumberOverlaps"] >= 2).drop("NumberOverlaps")
    print(f"Number of union terminal exons with >=2 overlapping PAS - {m_t_exons.te_id_min.nunique()}")

    # Filter & annotate PAS overlapping TEs with APA
    print(f"Number of unique PAS in input BED - {pas.pas_id.nunique()}")
    # Default is inner join so all PAS not overlapping TEs are dropped
    pas = pas.join(m_t_exons, strandedness="same").drop(["Start_b", "End_b", "Strand_b"])
    print(f"Number of unique PAS overlapping TEs with APA - {pas.pas_id.nunique()}")

    # Compute sum of TPMs of overlapping PAS for each terminal exon
    pas = pas.as_df()

    m_t_exons_sum = group_sum(pas,
                              group_col="te_id_min",
                              value_col="tpm",
                              out_col="sum_tpm")

    # print(m_t_exons_sum)

    # Select two representative sites for each TE by largest TPM
    # Note: Many TEs will only have 2 overlapping PAS so this filtering step is null. Opting to consider all events together for code simplicity
    pas_rep = pas.groupby("te_id_min").apply(lambda x: x.nlargest(2, "tpm")).reset_index(drop=True)

    # print(pas_rep)

    # Calculate sum of expression of two representative sites for each TE
    pas_rep_sum = group_sum(pas_rep,
                        group_col="te_id_min",
                        value_col="tpm",
                        out_col="sum_tpm")

    # print(pas_rep_sum)

    # Compute fraction of total expression that originates from two representative sites
    comb_sum = (pas_rep_sum.merge(m_t_exons_sum,
                                  on="te_id_min",
                                  suffixes=["_rep", "_tot"])
                .assign(rep_total_expr_frac=lambda df: df["sum_tpm_rep"] / df["sum_tpm_tot"])
                )

    # Get a set of terminal exons where representative fraction >= min_total_expr_frac of total expression on TE
    min_total_tes = set(comb_sum.loc[comb_sum["rep_total_expr_frac"] >= min_total_expr_frac,
                                     "te_id_min"])

    print(f"Number of union terminal exons where 2 representative sites represent at least {min_total_expr_frac} of total PAS expression on TE - {len(min_total_tes)}")

    # Subset representative PAS for TEs passing this filter
    pas_rep = pas_rep.loc[pas_rep["te_id_min"].isin(min_total_tes), :]

    # Check that 'minor PAS' for each TE has >= min_frac_site of total expression on TE (of representative sites)

    # First calculate relative expression (fraction) of each PAS on terminal exon
    pas_rep.loc[:, "rel_exp"] = (pas_rep.groupby("te_id_min",
                                                 group_keys=False)
                                 ["tpm"]
                                 .apply(lambda x: x / x.sum()))

    # print(pas_rep)

    # Filter for TEs where minor site has fractional expression >= min_frac_site
    pas_rep = pas_rep.groupby("te_id_min").filter(lambda x: x["rel_exp"].min() >= min_frac_site)
    print(f"Number of union terminal exons where minor site has >= {min_frac_site} of total expression on terminal exon - {pas_rep['te_id_min'].nunique()}")


    # Check that selected sites are not within window_size of one another (prevent matches to both prox & distal PAS)

    # First assign 'pas_number' to differentiate proximal and distal pas
    # 1 = proximal pas, 2 = distal pas
    pas_rep = pr.PyRanges(pas_rep, int64=True)
    # Have to assign a 'dummy column' to conform to add_region_number
    pas_rep = filterPAS.add_region_number(pas_rep.assign("Feature",
                                               lambda df: pd.Series(["pas"]*len(df.index))),
                                id_col="te_id_min",
                                out_col="pas_number",
                                feature_key="pas"
                                ).drop("Feature")

    # print(pas_rep[["te_id_min", "pas_number"]])

    # Overlap join proximal sites with distal sites
    # slack extends intervals prior to overlap query (i.e. allows overlap within window_size)
    # non-overlapping intervals have cols from query filled with -1
    pas_rep_prox = pas_rep.subset(lambda df: df["pas_number"] == 1)
    prox_dist = pas_rep_prox.join(pas_rep.subset(lambda df: df["pas_number"] == 2),
                                  how="left",
                                  strandedness="same",
                                  slack=window_size)

    # Subset to overlapping PAS
    prox_dist = prox_dist.subset(lambda df: df["Start_b"] != -1)

    if len(prox_dist) != 0:
        # Some prox and distal PAS overlap - check they come from same te_id
        prox_dist = prox_dist.subset(lambda df: df["te_id_min"] == df["te_id_min_b"])
        olap_ids = set(prox_dist.te_id_min)

        print(f"Number of union terminal exons where representative PAS fall within {window_size} of one another - {len(olap_ids)}")

        # Remove overlapping TEs
        pas_rep = pas_rep.subset(lambda df: ~df["te_id_min"].isin(olap_ids))


    else:
        print(f"Number of union terminal exons where representative PAS fall within {window_size} of one another - 0")

    print(f"Total number of TEs with APA after filtering steps - {pas_rep.te_id_min.nunique()}")

    # Generate full terminal exon ID & output BEDs
    # <gene_id>|<tx_id1>,<tx_id2>|<prediction_group_identifier>|<pas_number>
    # PAS - <gene_id>|<tx_id1>,<tx_id2>|NA|<1/2>
    # TEs - <gene_id>|<tx_id1>,<tx_id2>|NA|NA

    valid_te_id_min = set(pas_rep.te_id_min)

    # Generate te_id for PAS BED
    pas_rep.predn_group = "NA"

    pas_rep = filterPAS.assign_id(pas_rep,
                                  cols_to_cat=["te_id_min", "predn_group", "pas_number"],
                                  sep_char="|",
                                  out_col="te_id")

    pas_rep = pas_rep[["te_id", "rel_exp"]].apply(lambda df: df.rename(columns={"te_id": "Name",
                                                                                "rel_exp": "Score"}))

    print(pas_rep)

    # te_id for terminal exons BED
    m_t_exons = m_t_exons.subset(lambda df: df["te_id_min"].isin(valid_te_id_min))

    # Create dummy columns
    m_t_exons.predn_group = "NA"
    m_t_exons.pas_number = "NA"

    m_t_exons = filterPAS.assign_id(m_t_exons,
                                    cols_to_cat=["te_id_min", "predn_group", "pas_number"],
                                    sep_char="|",
                                    out_col="te_id")

    m_t_exons = m_t_exons[["te_id"]].apply(lambda df: df.rename(columns={"te_id": "Name"}))

    # Output to BED
    # Note: pr.to_bed will fill default cols if not present, so no need to manually add Score for terminal exons
    pas_rep.to_bed(out_prefix + ".rep_gt_pas.bed")
    m_t_exons.to_bed(out_prefix + ".rep_tes.bed")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]),sys.argv[6])
