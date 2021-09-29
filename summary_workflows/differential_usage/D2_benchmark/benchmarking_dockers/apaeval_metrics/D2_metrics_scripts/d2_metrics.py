#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 12:00:14 2021

@author: weronika
"""

import pandas as pd

def label_match(row_pd, row_gt):
    """
    Returns row labels depending on predicted and GT values.
    """

    if row_pd == 1 and row_gt == 1:
        return "TP"
    if row_pd == 1 and row_gt == 0:
        return "FP"
    if row_pd == 0 and row_gt == 0:
        return "TN"
    if row_pd == 0 and row_gt == 1:
        return "FN"
    else:
        return "None"
    
def count_matches(merged_df):
    """
    Returns counts of true positives, true negatives, false positives and false negatives.
    """
    
    matches = list(merged_df["match"])
    return matches.count("TP"), matches.count("TN"), matches.count("FP"), matches.count("FN"), 

def sensitivity(tp, fn):
    """
    Returns sensitivity of prediction.
    """
    
    return tp / (tp + fn)

def fdr(tp, fp):
    """
    Returns False Discovery Rate of prediction.
    """
    
    return fp / (tp + fp)  
    
    
def calculate_metrics(f_PD, f_GT):
    """
    Returns Sensitivity and False Discovery Rate of prediction.
    """
    
    df_PD = pd.read_csv(f_PD, header=None)
    df_PD.columns = ["gene", "is_de_PD", "is_lengthened_PD"]
    df_GT = pd.read_csv(f_GT, header=None)
    df_GT.columns = ["gene", "is_de_GT", "is_lengthened_GT"]

    merged_df = pd.merge(df_PD, df_GT, 'outer')
    # if gene data not available in either PD or GT, values set to "no differential expression detected":
    merged_df.fillna({'is_de_PD':0, 'is_lengthened_PD':0, "is_de_GT":0, "is_lengthened_GT":0}, inplace=True)
    merged_df = merged_df.astype({'is_de_PD':int, 'is_lengthened_PD':int, "is_de_GT":int, "is_lengthened_GT":int})
  
    merged_df["match"] = merged_df.apply(lambda row: label_match(row.is_de_PD, row.is_de_GT), axis=1)
   
    tp, tn, fp, fn = count_matches(merged_df)
    
    return sensitivity(tp, fn), fdr(tp, fp)
    
    
    
    
    
    
    
    
