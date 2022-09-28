import subprocess
import os
from io import StringIO
import pandas as pd
import pyranges as pr
import numpy as np
from scipy.stats import pearsonr, spearmanr
from shutil import which

#############################
# Genome/annotation related
#############################

def select_genome_file(file_name, genome_path):
    """Select the genome file according to the organism.

    Requires that the file_name contains an expression containing organism
    information, which will be matched against the genome_path directory.
    The format should be: name.mm10.ext or name.hg38extension.ext, with
    matching genome annotations: gencode.mm10.gtf and gencode.hg38extension.gtf.
    Note: no check for the extension (e.g. gtf) is done.

    Args:
        file_name (str): Name containing organism information. Supported: mm* and hg*.
        genome_path (str): directory containing genome annotations in gtf format.

    Returns:
        str with genome file path.
    """
    GENOME_STRINGS = ["mm", "hg"]
    SPLITSTRING = "."
    assert os.path.exists(genome_path), f"Genome annotation directory not found: {genome_path}"
    file_components =  file_name.split(SPLITSTRING)
    # search for genome
    for genome_string in GENOME_STRINGS:
        match = [comp for comp in file_components if genome_string in comp]
        if len(match) != 0:
            break
    if len(match) == 0:
        raise ValueError(f"No genome string: {GENOME_STRINGS} in file_name: {file_name} found.")
    # find all genome files in genome_path
    for f in os.listdir(genome_path):
        # find exact match in file
        genome_match = [f for comp in f.split(SPLITSTRING) if match[0] == comp]
        if len(genome_match) != 0:
            break
    if len(genome_match) == 0:
        raise ValueError(f"No genome string: {GENOME_STRINGS} in genome_path: {genome_path} found.")
    # return path to file
    return os.path.join(genome_path, genome_match[0])


def load_genome(genome_path):
    """Load genome annotation in gtf format.
    
    Requires feature 'gene', which is available in gencode.
    
    Args:
        genome_path: file path for gtf.
    
    Returns:
        pandas.df with genes.
    """
    assert os.path.exists(genome_path), f"Genome annotation not found: {genome_path}"
    # load gtf according to specifications (https://www.ensembl.org/info/website/upload/gff.html)
    genome = pd.read_csv(genome_path, sep="\t", header=None,
            comment="#",
            names=['seqname','source', 'feature','start','end','score','strand','frame','attribute'],
            usecols=['seqname','start','end','strand','feature'],
            dtype={'seqname': str, 'start': np.int64, 'end': np.int64, 'strand': str, 'feature': str})
    # Subset genome to only features gene
    subset = genome['feature'] == "gene"
    assert not subset.sum() == 0, f"No feature: 'gene' found in {genome_path}"
    genome.drop('feature', 1, inplace=True)
    return genome[subset]
    

def find_pas_genes(gene, df):
    """Find all PAS in given gene.

    Args:
        gene (dict): dict with keys 'seqname', 'start', 'end' and 'strand'.
        df (pandas.df): df with columns 'chrom_g' (str), 'chromStart_g' (int),
            'chromEnd_g' (int) and 'strand_g'.

    Returns:
        pd.Series: Array of indices indicating PAS for gene 'gene'.
    """
    subset = (df['chrom_g'] == gene['seqname']) & (df['strand_g'] == gene['strand'])
    subset = subset & (df['chromStart_g'] > gene['start']) & (df['chromEnd_g'] < gene['end'])
    return subset

#####################
# PAS filtering
#####################
## filtering PAS from ground truth
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

## filtering PAS from the prediction


#####################
# PAS matching
#####################


def bedtools_window(bed1, bed2, window, reverse=False):
    """
    Python wrapper for bedtools window.
    reverse: return only sites that have no match in the ground truth.
    """
    
    # make sure bedtools can be called in current env
    assert which('bedtools') is not None, "bedtools not installed or not in PATH"

    # run bedtools window, capture output
    if not reverse:
        out = subprocess.run(['bedtools', 'window', '-sm',
                              '-w', str(window), 
                              '-a', bed1, 
                              '-b', bed2], 
                             capture_output=True, shell=False)
    else:
        out = subprocess.run(['bedtools', 'window', '-sm', '-v',
                              '-w', str(window), 
                              '-a', bed1,
                              '-b', bed2], 
                             capture_output=True, shell=False)

    assert out.returncode==0, "bedtools window run failed, check input files"

    # memory file-handle to pass output to pandas without writing to disk
    out_handle = StringIO(out.stdout.decode())
    
    # incase there were no sites returned (no overlap / all overlap in case of reverse=True)
    if not out.stdout.decode():
        out = pd.DataFrame()
    else:
        out = pd.read_csv(out_handle, delimiter='\t', header=None, dtype={0: str,6: str})
        
    # label columns
    out.rename({0: 'chrom_p', 1: 'chromStart_p', 2: 'chromEnd_p', 3: 'name_p', 4: 'score_p', 5: 'strand_p', 6: 'chrom_g', 7: 'chromStart_g', 8: 'chromEnd_g', 9: 'name_g', 10: 'score_g', 11: 'strand_g'}, axis=1, inplace=True)
   
    return(out)

def find_weights(matched_sites_orig, pd_pos):
    """
    Calculate score weights for PD sites matched to multiple GT sites.
    If PD and GT sites are a perfect match, the weights are not split.
    If PD site between two GT sites, weights calculated based on distance.
    """
    matched_sites = matched_sites_orig.copy()
    matched_sites.sort_values(by=["dist"], inplace=True)
    matched_sites = matched_sites.reset_index(drop=True)
    matched_sites["weight"] = 0.0

    # if the PD is not between GT sites, assume the closest GT is a perfect match
    if pd_pos < matched_sites["chromEnd_g"].min() or pd_pos > matched_sites["chromEnd_g"].max():
        matched_sites["dist"][matched_sites["dist"].idxmin()] = 0.0

    # calculate weight based on distance; a perfect match will have weight set to 1.0
    matched_sites["weight"][0:2] = 1 - (matched_sites["dist"][0:2] / np.sum(matched_sites["dist"][0:2])) # weight is equal to distance proportion
            
    return matched_sites

def split_pd_by_dist(matched_sites):
    """
    Identify PD sites matched to multiple GT sites and split the score between 2 sites based on distance
    """
    # find rows in sample that are not unique 
    # -> extended predicted site overlaps multiple ground truth sites
    
    unique_pd_cols = ['chrom_p', 'chromStart_p', 'chromEnd_p', 'strand_p']
    not_unq_mask = matched_sites.duplicated(unique_pd_cols, keep=False)
    

    if np.sum(not_unq_mask) > 0: # otherwise there was no overlap
            
        not_unq = matched_sites[not_unq_mask].copy()
       
        # remove non-unique sites
        matched_sites.drop(matched_sites.index[not_unq_mask], inplace=True, axis=0)
        
        # calculate distances between PD and GT sites
        not_unq["dist"] = abs(not_unq["chromEnd_p"] - not_unq["chromEnd_g"])

         # calculate weights based on the distance to two nearest GT sites
        not_unq = not_unq.groupby(unique_pd_cols).apply(lambda x: find_weights(x, x.name[2])).reset_index(drop=True)
        # apply weights to PD score
        not_unq["score_p"] = not_unq["score_p"] * not_unq["weight"]
        # rename sites that have been split
        not_unq["name_p"] = not_unq["name_p"] + "_" + not_unq["name_g"]

        # merge back to main df
        not_unq.drop(["dist", "weight"], axis=1, inplace=True)
        matched_sites = pd.concat([matched_sites, not_unq])

    return matched_sites

def merge_pd_by_gt(matched_sites):
    """
    Identify multiple PD sites matched to single GT site, merge and sum the scores.
    """
    
    # find GT rows in sample that are not unique - multiple PD matched to one GT
    # -> extended predicted site overlaps multiple ground truth sites
    unique_gt_cols = ['chrom_g', 'chromStart_g', 'chromEnd_g', 'strand_g']
    not_unq_GT_mask = matched_sites.duplicated(unique_gt_cols, keep=False)

    if np.sum(not_unq_GT_mask) > 0: # otherwise there was no overlap
        
        not_unq_GT = matched_sites[not_unq_GT_mask].copy()
        # remove non-unique sites
        matched_sites.drop(matched_sites.index[not_unq_GT_mask], inplace=True, axis=0)
        
        # sum scores for all PD sites matched to the same GT
        not_unq_GT["score_p"] = not_unq_GT.groupby(unique_gt_cols)["score_p"].transform("sum")
        # use only GT location for merged PD sites
        not_unq_GT[["chrom_p", "chromStart_p", "chromEnd_p"]] = not_unq_GT[["chrom_g", "chromStart_g", "chromEnd_g"]]
        not_unq_GT["name_p"] = not_unq_GT["name_g"] + "_merged"

        # drop duplicates and merge back to main df
        not_unq_GT.drop_duplicates(subset=["name_g"], keep="first", inplace=True)
        matched_sites = pd.concat([matched_sites, not_unq_GT])

    return matched_sites

###############
# Metrics
###############


def precision(tp, fp):
    """
    Returns precision of prediction.
    """
    return tp / (tp + fp)

def sensitivity(tp, fn):
    """
    Returns sensitivity/TPR/recall of prediction.
    """
    return tp / (tp + fn)

def fdr(tp, fp):
    """
    Returns False Discovery Rate of prediction.
    """
    return fp / (tp + fp) 

def f1_score(precision, sensitivity):
    """
    Returns F1 score of prediction.
    """
    return  2 * (precision * sensitivity) / (precision + sensitivity)

def jaccard(tp, fp, fn):
    """
    Returns Jaccard index of prediction.
    """
    return  tp / (tp + fp + fn)

def corr_Pearson_with_gt(matched_sites):

    vec_true = matched_sites["score_g"]
    vec_pred = matched_sites["score_p"]

    # correlation coefficient
    r = pearsonr(vec_true, vec_pred)[0]

    return(r)

def corr_Spearman_with_gt(matched_sites):

    vec_true = matched_sites["score_g"]
    vec_pred = matched_sites["score_p"]

    # correlation coefficient
    r = spearmanr(vec_true, vec_pred)[0]

    return(r)

def relative_pas_usage(merged_bed_df, genome):
    """Compute correlation of relative PAS usage.

    1. find all PAS for a given gene, i.e. all ground truth PAS for a given gene or
        all predicted PAS that are matched to ground truth for a given gene.
    2. sum TPM values of all PAS
    3. calculate fraction for each PAS by dividing TPM_PAS by TPM_sum.

    Args:
        merged_bed_df (pandas.df): Table of matched prediction and ground truth PAS. From match_with_gt().
        genome (pandas.df): Table with gene positions.

    Returns:
        float: df of normalised and relative PAS values.
    """
    df = merged_bed_df.copy()
    # get list of PAS per gene
    pas_per_gene = [find_pas_genes(gene, df) for i, gene in genome.iterrows()]
    # remove PAS not covered by gene
    pas_per_gene = [pas for pas in pas_per_gene if len(np.where(pas)[0]) != 0]
    # compute fraction of PAS per gene
    normalised_dfs = [fraction_pas(gene_pas, df) for gene_pas in pas_per_gene]
    # concatenate list of pandas.df
    normalised_df = pd.concat(normalised_dfs, axis=0)
    return normalised_df

def fraction_pas(gene_pas, df):
    """Compute fraction for each PAS.
    Each PAS in the given gene is normalised by the sum of TPM values (column 'score'),
    separately for prediction and ground truth.

    Args:
        gene_pas (pandas.Series): bool vector indicating PAS for given gene.
        df (pandas.df): Table of matched prediction and ground truth PAS. From match_with_gt().
                        Columns 'score_p' and 'score_g' are used.
    
    Returns:
        pandas.df: 'df' for rows 'gene_pas' with 'score' columns as fractions.
    """
    pred_sum = df.loc[gene_pas, 'score_p'].sum()
    if pred_sum == 0:
        df.loc[gene_pas,'score_p'] = 0
    else:
        df.loc[gene_pas,'score_p'] = df.loc[gene_pas, 'score_p'] / pred_sum
    ground_sum = df.loc[gene_pas, 'score_g'].sum()
    if ground_sum == 0:
        df.loc[gene_pas,'score_g'] = 0
    else:
        df.loc[gene_pas,'score_g'] = df.loc[gene_pas,'score_g'] / ground_sum
    return df.loc[gene_pas,:]

