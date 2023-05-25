Module main
===========

Functions
---------

    
`bedtools_window(bed1, bed2, window, reverse=False)`
:   Python wrapper for bedtools window.
    reverse: return only sites that have no match in the ground truth.

    
`corr_Pearson_with_gt(matched_sites)`
:   

    
`corr_Spearman_with_gt(matched_sites)`
:   

    
`f1_score(precision, sensitivity)`
:   Returns F1 score of prediction.

    
`fdr(tp, fp)`
:   Returns False Discovery Rate of prediction.

    
`find_pas_genes(gene, df)`
:   Find all PAS in given gene.
    
    Args:
        gene (dict): dict with keys 'seqname', 'start', 'end' and 'strand'.
        df (pandas.df): df with columns 'chrom_g' (str), 'chromStart_g' (int),
            'chromEnd_g' (int) and 'strand_g'.
    
    Returns:
        pd.Series: Array of indices indicating PAS for gene 'gene'.

    
`find_weights(matched_sites_orig, pd_pos)`
:   Calculate score weights for PD sites matched to multiple GT sites.
    If PD and GT sites are a perfect match, the weights are not split.
    If PD site between two GT sites, weights calculated based on distance.

    
`fraction_pas(gene_pas, df)`
:   Compute fraction for each PAS.
    Each PAS in the given gene is normalised by the sum of TPM values (column 'score'),
    separately for prediction and ground truth.
    
    Args:
        gene_pas (pandas.Series): bool vector indicating PAS for given gene.
        df (pandas.df): Table of matched prediction and ground truth PAS. From match_with_gt().
                        Columns 'score_p' and 'score_g' are used.
    
    Returns:
        pandas.df: 'df' for rows 'gene_pas' with 'score' columns as fractions.

    
`jaccard(tp, fp, fn)`
:   Returns Jaccard index of prediction.

    
`load_genome(genome_path, feature='gene')`
:   Load genome annotation in gtf format.
    
    And subset to only return feature.
    E.g. feature 'gene', which is available in gencode.
    
    Note:
        The returned dataframe contains all attributes as read in by pyranges.
    
    Args:
        genome_path: file path for gtf.
        feature (str): feature to subset. Default: gene.
    
    Returns:
        pandas.df with [feature] rows from genome.

    
`merge_pd_by_gt(matched_sites)`
:   Identify multiple PD sites matched to single GT site, merge and sum the scores.

    
`precision(tp, fp)`
:   Returns precision of prediction.

    
`relative_pas_usage(merged_bed_df, genome)`
:   Compute correlation of relative PAS usage.
    
    1. find all PAS for a given gene, i.e. all ground truth PAS for a given gene or
        all predicted PAS that are matched to ground truth for a given gene.
    2. sum TPM values of all PAS
    3. calculate fraction for each PAS by dividing TPM_PAS by TPM_sum.
    
    Args:
        merged_bed_df (pandas.df): Table of matched prediction and ground truth PAS. From match_with_gt().
        genome (pandas.df): Table with gene positions.
    
    Returns:
        float: df of normalised and relative PAS values.

    
`select_genome_file(file_name, genome_path)`
:   Select the genome file according to the organism.
    
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

    
`sensitivity(tp, fn)`
:   Returns sensitivity/TPR/recall of prediction.

    
`split_pd_by_dist(matched_sites)`
:   Identify PD sites matched to multiple GT sites and split the score between 2 sites based on distance