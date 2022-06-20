#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 12:13:58 2021

@author: weronika
"""

import d2_metrics as metrics

file1 = "example_files/input1.csv"
file2 = "example_files/input2.csv"

s, f = metrics.calculate_metrics(file1, file2)

print(s, f)