# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:44:10 2016

@author: shooper
"""
import pandas as pd
import numpy as np

txt = '/vol/v2/stem/ebird/samples/z_leucophrys_res1yr_20161201_1550/erd_v5_Zonotrichia_leucophrys_predictors.txt'
df = pd.read_csv(txt, sep='\t', index_col='obs_id')
unique = df[['row','col']].drop_duplicates().values
stdv_obs = {}
stdv_time = {}
means = {}
for yr in range(2002, 2013):
    stdv_yr = []
    #stdv_yr_t = []
    #mean_yr = []
    df_yr = df[df.YEAR == yr]
    #unique = [(rc[0],rc[1]) for rc in df_yr[['row','col']].drop_duplicates().values]
    for r, c in unique:
        this_location = df_yr.ix[(df_yr.row == r) & (df_yr.col == c), 'z_leucophrys']
        if len(this_location) == 1:
            stdv = 0
            avg = this_location
            #stdv_t = 0
        else:
            stdv = this_location.std()
            avg = this_location.mean()
            #stdv_t = df_yr.ix[(df_yr.row == r) & (df_yr.col == c), 'TIME'].std()
        stdv_yr.append(stdv)
        #mean_yr.append(avg)
        #stdv_yr_t.append(stdv_t)
    stdv_obs[yr] = stdv_yr
    #stdv_time[yr] = mean_yr
    #stdv_time[yr] = stdv_yr_t

