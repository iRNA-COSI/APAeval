#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import sys

def get_terminal_exons(gr,
                       feature_col="Feature",
                       feature_key="exon",
                       id_col="transcript_id",
                       region_number_col="exon_number",
                       filter_single=False
                       ):
    '''
    Extract terminal exons for every transcript

    gr: pr.PyRanges object containing exons
    feature_col: Name of column in gr containing identifier of feature
    feature_key: feature identifier in feature_col. This must be the only value present in this column (i.e. only 'exons' are present)
    id_col: name of column by which to group regions (i.e. exons) in gr.
    region_number_col: name of column containing numbering regions (exons) of a group (transcript) in a strand-aware, 5'-3' ascending order
        - i.e. 1st exon (most 5') is numbered 1, terminal exon of transcript with n exons is labelled n
        - For reference GTFs this is generally a fair assumption
        - If not confident can use add_region_number to generate ranking fitting this expectation
        - must be an int np.dtype ('int32' or 'int64')
    filter_single: Whether to remove single exon transcripts (True) or not (False)
    '''

    assert region_number_col in gr.columns.tolist()
    assert feature_col in gr.columns.tolist()
    assert id_col in gr.columns.tolist()

    # Make sure region_number is an int dtype (so can be sorted numerically)
    assert gr.dtypes.loc[region_number_col] in [np.dtype('int32'), np.dtype('int64')], f"region_number_col - {region_number_col} - must be an int dtype"
    # Make sure only 'exon' features are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)


    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    mod_gr = gr.apply(lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True))


    # Filter out single-exon transcripts
    if filter_single:
        print("Filtering for multi-exon transcripts...")
        print(f"Before: {len(set(mod_gr.as_df()[id_col]))}")

        # Setting to 'False' marks all duplicates as True (so keeps transcript IDs with multiple exons these)
        mod_gr = mod_gr.subset(lambda df: df.duplicated(subset=[id_col], keep=False))

        print(f"After: {len(set(mod_gr.as_df()[id_col]))}")

    # Pick last region entry by max region number for each transcript (id_col)
    # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)

    out_gr = mod_gr.subset(lambda df: ~(df.duplicated(subset=[id_col], keep="last")))

    # re-sort by genomic coords not transcript_id (maybe unnecessary?)
    return out_gr.sort()


def _df_add_region_number(df, id_col, sort_col="Start"):
    '''
    Return a Series of strand-aware region numbers (5'-3' in 1..n)
    Function to be used internally in a pr.assign (mainly by add_region_number)
    '''
    if id_col not in df.columns.tolist():
        raise KeyError(f"id_col - {id_col} - is not present in df for chr/strand pair {','.join([df.Chromosome.iloc[0], df.Strand.iloc[0]])}")

    elif (df.Strand == "+").all():
        # Start position smallest to largest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=True)

    elif (df.Strand == "-").all():
        # Start position largest to smallest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=False)

    elif df.empty:
        print("df is empty - returning empty pd.Series")
        return pd.Series()


def add_region_number(gr,
                      id_col="transcript_id",
                      out_col="exon_number",
                      feature_col="Feature",
                      feature_key="exon",
                      ):
    '''
    Adds column to gr containing a strand aware region number column, ordered 5'-3' 1..n by a group of features (e.g. transcript)

    - Assumes every region/row of feature in group is not duplicated (i.e. given exon only appears once)
    - Assumes regions are not overlapping
        - Uses 'Start' coordinates to order regions from left-right
        - Overlapping regions with different 'Start' coordinates would be ranked different ranks
    '''

    # Make sure only 'feature_key' rows are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)
    assert gr.stranded

    # Make sure sorted by position first.
    gr = gr.sort()

    # Add in region number column in strand aware manner, so 1 = most 5', n = most 3'

    gr = gr.assign(out_col, lambda df: _df_add_region_number(df, id_col))

    return gr



def _collapse_tes(df,
                  cluster_col,
                  distinct_cols,
                  sep=","
                  ):
    '''
    Collapse grouped intervals (rows) into a single union interval whilst collapsing metadata
    Intended for internal pd.DataFrames of pr.PyRanges objects accessed through a pr.apply() call
    '''

    # Track input column order so can retain in final output
    col_order = df.columns.tolist()

    other_cols = [col for col in col_order if col not in distinct_cols]
    # First define how want to collapse metadata columns to single row values
    #
    # Cols want to combine values as ; separated string (e.g. if expect diff values per row)
    cat_dict = {col: lambda x: x.astype(str).str.cat(sep=sep)
                  for col in df.columns
                    if col not in distinct_cols}

    # Cols want to report a single value (e.g. where expect same values in all rows/not important)
    # Also important where need to maintain original
    distinct_dict = {col: lambda x: x.first()
                         for col in distinct_cols
                        }

    # Define the merged interval
    cat_dict["Start"] = lambda x: x.min()
    cat_dict["End"] = lambda x: x.max()

    grpd = df.groupby(cluster_col, group_keys=False)

    # Cols want to report a single value (e.g. where expect same values in all rows/ output value not important)
    grpd_distinct = grpd[distinct_cols].first()

    grpd_other = grpd[other_cols].agg(cat_dict)

    # Combine to single row
    clpsd = grpd_distinct.merge(grpd_other, left_index=True, right_index=True)

    return clpsd[col_order]


def collapse_tes(gr,
                 cluster_col="Cluster",
                 distinct_cols=["Feature", "gene_id"],
                 sep=","):
    '''
    Collapse overlapping regions into union interval whilst preserving metadata
    gr: pr.PyRanges object
    cluster_col: Name of column grouping intervals wish to collapse. e.g. output of pr.cluster()
    distinct_cols: list of column names wish to retain a single value for collapsed interval.
        - All other columns in object will be collapsed to a single string
        - Note: cluster_col, 'Chromosome' and 'Strand' will be added internally (required to produce a valid object)
    sep: string character to separate values of all columns not found in distinct_cols
    '''

    # Want to keep grouping col to single value
    # Chromosome and Strand must be a single value to output a valid PyRanges dataframe
    dist_cols = distinct_cols + [cluster_col] + ["Chromosome", "Strand"]

    return gr.apply(lambda df: _collapse_tes(df, cluster_col, dist_cols, sep=sep))


def assign_id(gr, cols_to_cat=["Chromosome","Start","End","Strand"], sep_char=":", out_col="pas_id"):
    '''
    Generate a identifier column from a subset of column names
    gr: pr.PyRanges object
    cols_to_cat: list of columns names wish to combine (order of names defines order in output column)
    sep_char: str to denote separator between values from columns defined in cols_to_cat
    out_col: str denoting name of output column
    '''

    gr_cols = gr.columns.tolist()

    for col in cols_to_cat:
        assert col in gr_cols, f"Following column name - {col} - not found in input gr"

    return gr.assign(out_col,
                     lambda df: df[cols_to_cat[0]].str.cat(df[cols_to_cat[1:]].astype(str),
                                                           sep=sep_char
                                                           )
                    )

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


def filterPAS(gt_bed,gtf,min_total_expr_frac, min_frac_site, window_size, out_prefix):
    '''
    '''

    # Read in GTF, subset to exons
    exons = pr.read_gtf(gtf).subset(lambda df: df["Feature"] == "exon")

    # extract terminal exons for each transcript
    exons = exons.assign("exon_number",
                         lambda df: df["exon_number"].astype(float).astype(int))

    t_exons = get_terminal_exons(exons)

    # Merge overlapping terminal exons into union intervals whilst retaining metadata

    # Subset to minimal columns
    t_exons = t_exons[["gene_id", "transcript_id"]]

    # print(t_exons)

    # First annotate overlapping intervals (on the same strand) with a common identifier
    # Only terminal exons of the same gene are grouped together in case of overlaps
    # Note: consider requiremnt for common gene_id. Any unintended consequences?
    t_exons = t_exons.cluster(strand=True, by="gene_id")

    # Collapse to union intervals, collapsing contributing transcript IDs to comma-separated string
    m_t_exons = collapse_tes(t_exons,
                                       distinct_cols=["gene_id"],
                                       sep=",")

    m_t_exons = m_t_exons.drop("Cluster")

    # print(m_t_exons)

    # Construct the minimal terminal exon ID - gene_id|<transcript_id1,transcript_id2>

    m_t_exons = assign_id(m_t_exons,
                                    cols_to_cat=["gene_id", "transcript_id"],
                                    sep_char="|",
                                    out_col="te_id_min")

    m_t_exons = m_t_exons.drop(["gene_id", "transcript_id"])

    # print(m_t_exons)

    # read in ground truth polyA sites
    pas = pr.read_bed(gt_bed)

    # assign an id column (pas_id) so can track number of PAS
    # Rename score column so clear it corresponds to TPM
    pas = assign_id(pas)
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
    pas_rep = add_region_number(pas_rep.assign("Feature",
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

    pas_rep = assign_id(pas_rep,
                                  cols_to_cat=["te_id_min", "predn_group", "pas_number"],
                                  sep_char="|",
                                  out_col="te_id")

    pas_rep = pas_rep[["te_id", "rel_exp"]].apply(lambda df: df.rename(columns={"te_id": "Name",
                                                                                "rel_exp": "Score"}))

    # print(pas_rep)

    # te_id for terminal exons BED
    m_t_exons = m_t_exons.subset(lambda df: df["te_id_min"].isin(valid_te_id_min))

    # Create dummy columns
    m_t_exons.predn_group = "NA"
    m_t_exons.pas_number = "NA"

    m_t_exons = assign_id(m_t_exons,
                                    cols_to_cat=["te_id_min", "predn_group", "pas_number"],
                                    sep_char="|",
                                    out_col="te_id")

    m_t_exons = m_t_exons[["te_id"]].apply(lambda df: df.rename(columns={"te_id": "Name"}))

    # Output to BED
    # Note: pr.to_bed will fill default cols if not present, so no need to manually add Score for terminal exons
    pas_rep.to_bed(out_prefix + ".rep_gt_pas.bed")
    m_t_exons.to_bed(out_prefix + ".rep_tes.bed")


if __name__ == '__main__':

    descrpn="""Filter a BED file of polyA sites to 2 representative sites overlapping terminal exons"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-b", "--bed",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to BED file of polyA sites")

    parser.add_argument("-g",
                        "--gtf",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to GTF file of transcript models")

    parser.add_argument("--min-total-expression-frac",
                        type=float,
                        default=0.8,
                        help="Minimum fraction of total expression of overlapping PAS provided by two highest expressed sites for terminal exon to be retained")

    parser.add_argument("--site-min-expression-frac",
                        type=float,
                        default=0.05,
                        help="Minimum fraction of total expression of representative sites on terminal exon contributed by minor PAS for terminal exon to be retained")

    parser.add_argument("-w", "--window-size",
                        required=True,
                        type=int,
                        default=50,
                        help="Size of window (nt) to consider polyA sites separated by less than window as overlapping")

    parser.add_argument("-o","--output-prefix",
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Name of prefix for output BED files. polyA site BED suffixed with '.rep_gt_pas.bed' & terminal exon BED '.rep_tes.bed'")


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # check that fractions are between 0-1
    assert 0 <= args.min_total_expression_frac <= 1, f"--min-total-expression-frac must be between 0 & 1, provided value - {args.min_total_expression_frac}"
    assert 0 <= args.site_min_expression_frac <= 1, f"--site-min-expression-frac must be between 0 & 1, provided value - {args.site_min_expression_frac}"

    main(args.bed,
         args.gtf,
         args.min_total_expression_frac,
         args.site_min_expression_frac,
         args.window_size,
         args.output_prefix)
