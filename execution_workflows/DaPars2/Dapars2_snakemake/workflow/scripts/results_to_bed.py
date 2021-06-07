#!/usr/bin/env python3


import pandas as pd
import numpy as np
import argparse
import os


def proximal_site_to_bed_coords(df,
                                coord_col="Predicted_Proximal_APA",
                                one_based=True,
                                start_outcol="start",
                                end_outcol="end"):
    '''
    Take the predicted proximal poly(A) site in column 'coord_col' (default: "Predicted_Proximal_APA")
    and return df with columns containing BED compliant start (named start_outcol) & end (named end_outcol) coordinates
    Coordinate is assumed to be one-based by default (one_based=True)

    Assuming the 'A' was the poly(A) site

    GFF      1 ┌─2─┐ 3   4   ...   - Start and end included (first base = 1)
             G---A---T   C   ...
    BED      0 └─1 └─2   3   ...   - Start included, end excluded (first base = 0)

    (representation stolen from https://gist.github.com/ilevantis/6d6ecf8718a5803acff736c2dffc933e#gff-vs-bed-indexing)
    '''

    return df.assign(**{start_outcol: lambda x: x[coord_col] - 1 if one_based
                        else x[coord_col],
                        end_outcol: lambda x: x[coord_col] if one_based
                        else x[coord_col] + 1
                        }
                    )


def distal_site_to_bed_coords(df,
                              coord_col = "Loci",
                              strand_col = "strand",
                              one_based = False,
                              start_outcol = "start",
                              end_outcol = "end"):
    '''
    Distal site will come from 3'end of provided annotation
    Since Dapars2 generates and uses BED files of annotated 3'UTRs,
    will by default assume starts are inclusive & ends exclusive (one_based = False)

    e.g. assuming a loci of (GFF - 1-2) and (BED - 0-2)


    GFF    ┌─1   2─┐ 3   4   ...
             G---A   T   C   ...
    BED    └─0   1 └─2   3   ...

    If + strand Then the poly(A) site (loci_end) becomes
    GFF      1 ┌─2─┐ 3   4   ...
             G   A   T   C   ...
    BED      0 └─1 └─2   3   ...

    if - strand then poly(A) site (loci_start) becomes
    GFF    ┌─1─┐ 2   3   4   ...
             G   A   T   C   ...
    BED    └─0 └─1   2   3   ...

    i.e. if going from 1-based (GFF) --> BED (one_based = True)
    + strand: start = loci_end - 1, end = loci_end
    - strand: start = loci_start - 1, end = loci_start

    if going from 0-based (BED) --> BED (one_based = False)
    + strand: start = loci_end - 1, end = loci_end (end coordinates in BED are exclusive)
    - strand: start = loci_start, end = loci_start + 1 (start coords in BED are inclusive)
    '''

    #First pull out start and end coords from 'Loci' column
    #formatted like chr6:119923969-119926687

    df[["loci_start", "loci_end"]] = (df[coord_col]
                                      .str.split(":", expand = True) # df with cols of chr6 & 119923969-119926687
                                      .loc[:,1] # 119923969-119926687 i.e. 2nd column
                                      .str.split("-", expand = True) #df with cols 119923969 & 119926687
                                      .astype(int)
                                     )

    # Select 3-prime most coordinate for each Tx from loci_start & end depending on strand
    df = df.assign(three_p = lambda x: np.where(x[strand_col] == "+",
                                                x["loci_end"],
                                                x["loci_start"]
                                               )
                  )

    # define BED start and end coords from 3-prime most coordinate
    if one_based:
        # GFF -> BED
        df = df.assign(**{start_outcol: lambda row: row["three_p"] - 1,
                          end_outcol: lambda row: row["three_p"]
                         }
                      )

    else:
        # BED of region to BED of single nucleotide
        df = df.assign(**{start_outcol: lambda x: np.where(x[strand_col] == "+",
                                                           x["three_p"] - 1,
                                                           x["three_p"]
                                                          ),
                          end_outcol: lambda x: np.where(x[strand_col] == "+",
                                                         x["three_p"],
                                                         x["three_p"] + 1 # - strand: 3'most = loci_start, which is inclusive in BED
                                                        )
                         }
                      )

    return df.drop(["loci_start","loci_end","three_p"], axis = 1)


def polya_site_id(df,
                  cols_to_combine = ["transcript_id", "gene_symbol"],
                  site_suffix = "prox",
                  sep = "_",
                  id_outcol = "name"):
    '''
    add an id column (named id_outcol),
    ID created by combining values in columns specified in cols_to_combine & a suffix (site_suffix)
    with values separated by '_'
    '''

    df["suffix"] = site_suffix

    cols = cols_to_combine + ["suffix"]

    df[id_outcol] = df[cols].apply(lambda row: sep.join(row.values.astype(str)), axis=1)

    return df.drop(["suffix"], axis = 1)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Script to generate BED files of predicted proximal and distal poly(A) sites from DaPars2 output tables")

    parser.add_argument("-i", "--input-tsvs", nargs="+", default=[], help = "space separated list of file paths to TSV files generated by DaPars2_Multi_Sample<_Multi_Chr>.py", required=True)
    parser.add_argument("-c","--coordinate-type", type=str, default="gff", choices=["bed","gff"], help="do regions specified in 'loci' column follow GFF or BED conventions? (default: %(default)s)")
    parser.add_argument("-o", "--output-bed", type=str, default="dapars2_identified_polya_sites.bed", help="name of output BED file (default: %(default)s)")

    args = parser.parse_args()

    flist = args.input_tsvs
    # print(flist)
    out_path = args.output_bed
    # print(out_path)

    if args.coordinate_type == "bed":
        loci_one_based = False
    elif args.coordinate_type == "gff":
        loci_one_based = True
    else:
        raise ValueError("'-c'/'--coordinate-type' must be one of 'bed' or 'gff' - {} was provided".format(args.coordinate_type))

    all_chrs = pd.concat([pd.read_csv(f, sep = "\t") for f in flist]).reset_index(drop=True)

    # 1. Pull out transcript_id, gene_name, chromosome and strand from region ID
    all_chrs[["transcript_id", "gene_symbol", "chromosome", "strand"]] = all_chrs["Gene"].str.split("|", expand = True)


    # 2. Pull out predicted proximal poly(A) site and convert to BED coordinates
    prox_all_chrs = proximal_site_to_bed_coords(all_chrs)

    # 3. Pull out distal site from the 'last exon' used for each transcript
    # e.g. loci column formatted like chr6:119923969-119926687
    distal_all_chrs = distal_site_to_bed_coords(all_chrs,one_based=loci_one_based)

    # 4. Assign an ID for each poly(A) site
    prox_all_chrs = polya_site_id(prox_all_chrs, site_suffix = "proximal")
    distal_all_chrs = polya_site_id(distal_all_chrs, site_suffix = "distal")

    # 5. BED requires a score column (5th) - assign empty field to '.' (identification BED)
    coords_all_chrs = pd.concat([prox_all_chrs, distal_all_chrs]).reset_index(drop=True)
    coords_all_chrs["score"] = "."

    # 6. Output to BED
    bed_cols = ["chromosome", "start", "end", "name", "score", "strand"]

    out_bed = coords_all_chrs[bed_cols]

    out_bed.to_csv(out_path, sep = "\t", header = False, index = False)
