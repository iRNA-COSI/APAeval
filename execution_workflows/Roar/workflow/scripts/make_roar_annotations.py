#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import sys

'''
Script to generate Roar-compliant GTF annotations
'''


# Pyranges default GTF columns are named as below
#https://github.com/biocore-ntnu/pyranges/blob/1ee215c645f7dbec3282555fcd0ccec610236614/pyranges/out.py#L44
pyranges_gtf_cols = "Chromosome   Source   Feature    Start     End       Score    Strand   Frame".split()

# Minimal cols needed from GTF for processing (no need to carry loads of cols)
processing_gtf_cols = pyranges_gtf_cols + ["gene_id", "transcript_id"]

# Cols required in output file
roar_gtf_outcols = pyranges_gtf_cols + ["gene_id", "length"]


def check_int64(gr1, gr2):
    '''
    Check Start and End columns in 2 PyRanges are 'int64' dtype and coerce if not
    '''

    gr1 = gr1.apply(lambda df: df.astype({"Start": "int64",
                                          "End": "int64"}
                                         )
                    )

    gr2 = gr2.apply(lambda df: df.astype({"Start": "int64",
                                          "End": "int64"
                                          }))

    return gr1, gr2


def longest_tx_per_gene(gr, gene_id_col="gene_id", tx_id_col="transcript_id"):
    '''
    '''

    # Get a df of gene_id | tx_id
    tx2gene = (gr.subset(lambda df: df["Feature"] == "transcript")
               .as_df()
               [[gene_id_col, tx_id_col]]
               .drop_duplicates())

    # Calculate transcript lengths by sum of lengths of individual exons
    # 'transcript' entries span introns + exons in genomic space
    exons = gr.subset(lambda df: df["Feature"] == "exon")
    exons.region_length = exons.lengths()

    # returns dict of {(chr, strand): df of <tx_id> | <tx_length>}
    tx_lengths = (exons.apply(lambda df: (df.groupby(tx_id_col)
                                          ["region_length"].sum()
                                          .reset_index()
                                          .rename({"region_length": "tx_length"},
                                                  axis=1)
                                          ),
                              as_pyranges=False
                              )
                  )

    # print(tx_lengths)
    tx_lengths = pd.concat(tx_lengths, ignore_index=True)

    # Get set of longest transcripts for each gene
    longest_txs = (set(tx2gene.merge(tx_lengths, on=tx_id_col)
                       .iloc[lambda df: df.groupby(gene_id_col)["tx_length"].idxmax(), :]
                       [tx_id_col]
                       )
                   )

    return gr.subset(lambda df: df[tx_id_col].isin(longest_txs))


def _sort_exons_by_tx(df, id_col):
    '''
    Intended to be applied inside a pr.apply() call (i.e. operates on underlying df of PyRanges object)
    '''

    assert df.Strand.nunique() == 1, f"df must contain only a single value in 'Strand' column, {df.Strand.nunique()} found"
    assert (df.Strand.drop_duplicates().isin(["+", "-"])).any(), f"'Strand' column values must be one of '+' or '-', {df.Strand.drop_duplicates().tolist()} found"

    if (df.Strand == "+").all():
        # Sort by transcript_id & exons 5' - 3' (smallest = left-most coord)
        # Ascending order = first in group is most 5' (first) exon
        return df.sort_values(by=[id_col, "Start"], ascending=True)

    elif (df.Strand == "-").all():
        # Sort by transcript_id & exons 5' - 3' (End coord = 5' end, largest = most 5')
        # Descending order = first in group is most 5' (first) exon
        return df.sort_values(by=[id_col, "Start"], ascending=False)


def annotate_last_exon(gr, id_col="transcript_id", out_col="last_exon"):
    '''
    Add a column specifying whether region corresponds to the most 3' region in the group (e.g. transcript) or not (1/0)
    '''

    # Sort grs by id_col & 5'-3' position (strand aware, first row of id is always most 5' region)
    gr = gr.apply(lambda df: _sort_exons_by_tx(df, id_col))

    # keep="last" sets last row by ID (last exon) to False and all others True
    gr = gr.assign(out_col,
                   lambda df: pd.Series(np.where(df.duplicated(subset=[id_col],
                                                               keep="last"),
                                                 0,
                                                 1)))

    return gr


def _most_tss_proximal(df, start_col="Start_b", end_col="End_b"):
    '''
    Return most 5' polyA site (strand-aware)
    '''
    assert start_col in df.columns
    assert end_col in df.columns
    assert df.Strand.nunique() == 1, f"df must contain only a single value in 'Strand' column, {df.Strand.nunique()} found"
    assert (df.Strand.drop_duplicates().isin(["+", "-"])).any(), f"'Strand' column values must be one of '+' or '-', {df.Strand.drop_duplicates().tolist()} found"

    if (df.Strand == "+").all():
        return df.loc[df[start_col].idxmin(), :]

    elif (df.Strand == "-").all():
        return df.loc[df[end_col].idxmax(), :]


def _representative_proximal_site(df, utr3_col="last_exon", utr3_key=1):
    '''
    intended for groupby pandas object
    '''

    utr_df = df.loc[df[utr3_col] == utr3_key, :]

    if len(utr_df) > 0:
        # Has 3'UTR overlapping, prioritise these
        return _most_tss_proximal(utr_df)

    else:
        # No 3'UTR overlapping, just select the most tss proximal of remaining
        return _most_tss_proximal(df)


def _df_change_coordinate(df, change_5p=False, suffix="_b"):
    '''
    Intended to be applied inside a pr.apply() call (i.e. operates on underlying df of PyRanges object)
    '''

    if (df.Strand == "+").all():
        if change_5p:
            # Start = 5'end
            # New Start = Start + suffix
            df["Start"] = df["Start" + suffix]
            return df

        else:
            # End = 3'end of exon
            # End + suffix = overlapping polyA site
            df["End"] = df["End" + suffix]
            return df

    elif (df.Strand == "-").all():
        if change_5p:
            # End = 5'end
            # New End = End + suffix
            df["End"] = df["End" + suffix]
            return df

        else:
            # Want to change 3' coordinate
            # Start = 3'end of exon
            # Start_b = Overlapping polyA site
            df["Start"] = df["Start" + suffix]
            return df


def get_pre_regions(gr, exons):
    '''

    Introns - stretch of DNA starting from the beginning of the exon containing
(or proximal to) the APA site and ending at the site position.
        - Need to find nearest upstream exon and use that 5'end
    '''

    gr_ex = gr.subset(lambda df: df.Feature == "exon")
    gr_int = gr.subset(lambda df: df.Feature == "intron")

    # Generate 'PRE' regions for exons (already have 5'end)
    # Update 3'end of feature to overlapping polyA site
    pre_ex = gr_ex.apply(lambda df: _df_change_coordinate(df,
                                                          change_5p=False,
                                                          suffix="_b"))

    # Generate 'PREP' regions for introns

    # Need to find nearest upstream exon to update 5'end
    # slack=1 reports bookended intervals as overlapping
    # NB: tried gr.nearest(how="upstream"), but nearest reported exon isn't always from the same gene

    # print(gr_int.columns)

    if len(gr_int) == 0:
        pre_int = pr.PyRanges()

    else:
        gr_int_up = gr_int.join(exons, slack=1, suffix="_c")

        # exons from other genes could be reported as overlapping
        gr_int_up = gr_int_up.subset(lambda df: df["gene_id"] == df["gene_id_c"])

        # slack=1 means downstream exon is also bookended (reported as overlapping)
        # Subset for the upstream exon (3'end = 5'end of intron)
        # (+ strand = End_c = Start)
        # (- strand = Start_c = End)
        gr_int_up = gr_int_up.subset(lambda df: (df.Strand == "+") & (df["Start"] == df["End_c"]) |
                                                (df.Strand == "-") & (df["End"] == df["Start_c"])
                                     )

        #     print(gr_int.gene_id.nunique())-
        #     print(gr_int_up.gene_id.nunique())
        #     print(gr_int_up.subset(lambda df: df.gene_id == df.gene_id_c).gene_id.nunique())

        # Now can generate PRE regions as before
        # First update the 5'end
        pre_int = gr_int_up.apply(lambda df: _df_change_coordinate(df, change_5p=True, suffix="_c"))
        # Now 3'end in same way as exons
        pre_int = pre_int.apply(lambda df: _df_change_coordinate(df, change_5p=False))

    return pr.concat([pre_ex, pre_int]).sort().drop(like="_c$")


def _pd_merge_gr(df, df_to_merge, how, on, suffixes, to_merge_cols):
    '''
    Perform a pd.merge inside a pr.apply to add columns from gr_to_merge based on metadata in gr_df
    Here, will check the chromosome and strand of provided gr (gr_df)
    and subset gr_to_merge to chr/strand before converting to df
    PyRanges dfs must have same columns across chr/strand pairs
    If not present in to merge need to output expected columns in df
    This should cut down on memory requirements vs convert to df (which requires pd.concat() x chrom & strand)
    For this kind of merge, only expect joins between regions on same chr/strand
    '''
    #chromsomes returns a list of chr names (always 1 val)
    assert isinstance(to_merge_cols, list)

    if df_to_merge.empty:
        print("df_to_merge for chr/strand pair {} is empty - returning to_merge cols filled with NaNs".format(",".join([df.Chromosome.iloc[0], df.Strand.iloc[0]])))

        df_cols = df.columns.tolist()
        # on won't have suffix added - need to remove as a target column
        to_merge_cols = [col for col in to_merge_cols if col != on]

        # print(df_cols)
        # print(to_merge_cols)

        # list of cols shared between dfs - need to add suffixes[1]
        # list of cols only in df_to_merge - these stay the same in a merge
        only_cols = [col for col in to_merge_cols if col not in df_cols]
        shared_cols = [col for col in to_merge_cols if col in df_cols]

        # print("only_cols - {}".format(only_cols))
        # print("shared_cols - {}".format(shared_cols))

        out_shared = [col + suffixes[1] for col in shared_cols]
        target_cols = out_shared + only_cols

        nrows = len(df.index)

        out_cols = {col: pd.Series([np.nan]*nrows) for col in target_cols}

        return df.assign(**out_cols)

    else:
        return df.merge(df_to_merge,
                        how=how,
                        on=on,
                        suffixes=suffixes)


def get_post_regions(pre_regions, tx_gr):
    '''
    '''

    #1. Get canonical transcript ends
    ends_3p = tx_gr.three_end()

    ends_cols = ends_3p.columns.tolist()

#     print(pre_regions.gene_id.nunique())

    #2. Join pre regions with 3'ends (coords on same row)
    pre_regions_3p = pre_regions.apply_pair(ends_3p,
                                            lambda pre_reg, ends_3p: _pd_merge_gr(pre_reg,
                                                                                  ends_3p,
                                                                                  how="inner",
                                                                                  on="gene_id",
                                                                                  suffixes=[None, "_3p"],
                                                                                  to_merge_cols=ends_cols))

#     print(pre_regions_3p.gene_id.nunique())
#     print(pre_regions)

    #3. Update Pre region End ('End') to canonical tx end
    pre_regions_3p = pre_regions_3p.apply(lambda df: _df_change_coordinate(df, change_5p=False, suffix="_3p"))

    #4. Update

    #4. Subtract PRE region from extended region to get POST region
    pre_regions_3p, pre_regions = check_int64(pre_regions_3p, pre_regions)

    # print(pre_regions.dtypes)
    # print(pre_regions_3p.dtypes)

    post_regions = pre_regions_3p.subtract(pre_regions)

    return post_regions.drop(like="_3p$")


def main(gtf_path,
         polya_bed_path,
         out_gtf_path):
    '''
    '''


    gtf = pr.read_gtf(gtf_path)
    pas = pr.read_bed(polya_bed_path)

    gtf = gtf[processing_gtf_cols]

    #1. Get set of gene IDs with at least 1 overlapping pas
    apa_gene_ids = set(gtf.overlap(pas).gene_id)

    gtf_apa = gtf.subset(lambda df: df["gene_id"].isin(apa_gene_ids))

    #2. For each gene, select the longest transcript as representative model
    gtf_apa_lng = longest_tx_per_gene(gtf_apa)

    # Generate separate 'features' objects for these transcripts
    exon_apa_lng = gtf_apa_lng.subset(lambda df: df["Feature"] == "exon")
    intron_apa_lng = gtf_apa_lng.features.introns(by="transcript")
    tx_apa_lng = gtf_apa_lng.subset(lambda df: df["Feature"] == "transcript")

    #3. Add column to annotate last exons (1/0)
    # Roar prioritises overlapping polyA sites in 3'UTRs
    # Gencode/Ensembl annotate UTR features differently
    # Focusing on exons is more general

    exon_apa_lng = annotate_last_exon(exon_apa_lng)

    #4. Select representative proximal polyA site for each tx
    # Prioritising overlaps in 3'UTRs (last exons here):
    # - pick overlapping pA site closest to transcription start site (most 5')

    prox_sites = (pr.concat([exon_apa_lng, intron_apa_lng])
                  # Add overlapping pA to same row (with _b suffix for cols)
                  .join(pas)
                  # Pick most 5' overlapping pA (End_b on plus, Start_b on minus)
                  .apply(lambda df: (df.groupby("transcript_id")
                                     .apply(lambda x: _representative_proximal_site(x))
                                     )
                         )
                  )

    #4. Define 'PRE' regions
    # Starts at 5'end of overlapping exon, ends at polyA site
    # If in intron, find upstream exon and use this exon as the 5'end
    pre_regions = get_pre_regions(prox_sites, exon_apa_lng)

    #5. Define 'POST' regions
    # Starts at 3'end of PRE region, ends at annotated tx end
    post_regions = get_post_regions(pre_regions, tx_apa_lng)

    #6. Roar-compliant gene IDs
    # PRE regions should be suffixed with '_PRE'
    # POST regions should be suffixed with '_POST'
    # Also add a 'length' attribute reporting the length (nt) of each region

    pre_regions.length = pre_regions.lengths()
    pre_regions = pre_regions.assign("gene_id",
                                     lambda df: df["gene_id"] + "_PRE")

    post_regions.length = post_regions.lengths()
    post_regions = post_regions.assign("gene_id",
                                       lambda df: df["gene_id"] + "_POST")

    #7. Combine into single PyRanges and output GTF
    all_regions = pr.concat([pre_regions[roar_gtf_outcols],
                             post_regions[roar_gtf_outcols]
                             ]).sort()

    all_regions.to_gtf(out_gtf_path)



if __name__ == '__main__':

    descripn = """Script to generate Roar compliant GTF annotations ('single APA' version) from reference GTF and BED file of polyA sites"""
    parser = argparse.ArgumentParser(description=descripn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-g",
                        "--gtf",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to reference GTF file (MUST contain 'gene_id' & 'transcript_id' attributes)")

    parser.add_argument("-p",
                        "--polya",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to BED file of reference polyA sites. 1nt length coordinates are strongly recommended")

    parser.add_argument("-o",
                        "--output-gtf",
                        default="roar_single_apa.gtf",
                        type=str,
                        help="Path to/name of output GTF file")

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    main(args.gtf, args.polya, args.output_gtf)
