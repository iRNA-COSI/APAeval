import pyranges as pr
import pandas as pd
import numpy as np


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
