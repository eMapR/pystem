# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 21:35:22 2016

@author: shooper
"""

import os
import sys
import shutil
import fnmatch
import pandas as pd
import numpy as np
from gdalconst import *
from scipy import stats
import time
import stem

def read_params(txt):
    '''
    Return a dictionary and a dataframe from parsed parameters in txt
    '''
    if not os.path.exists(txt):
        print 'Param file does not exist:\n', txt
    d = {}
    
    # Read in the rest of the text file line by line
    try:
        with open(txt) as f:
            input_vars = [line.split(";") for line in f]     
    except: 
        print 'Problem reading parameter file:\n', txt
        return None
    
    # Set the dictionary key to whatever is left of the ";" and the val
    #   to whatever is to the right. Strip whitespace too.
    n_skip_lines = 0 #Keep track of the number of lines w/ a ";"
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '{0}'.format(var[1].strip(" ").replace("\n", ""))
            n_skip_lines +=1
    
    # Get the lines with information about each variable as a df
    ''' Consider making skip_lines a list of indices rather than a range (maybe search \t)'''
    skip_lines = range(len(input_vars) - n_skip_lines, len(input_vars))
    df_vars = pd.read_csv(txt, sep='\t', index_col='var_name', skip_blank_lines=True, skiprows=skip_lines)
    # Drop any rows for which basepath or search str are empty
    df_vars.dropna(inplace=True, subset=['basepath','search_str'])
    df_vars.fillna({'path_filter': ''}, inplace=True)
    
    print '\nParameters read from:\n', txt, '\n'
    return d, df_vars


def get_predictors(years, search_dir, search_str, obs_txt, index_col, year_col, df_vars):
       
    df_obs = pd.read_csv(obs_txt, sep='\t', index_col=index_col)
    
    # Get a dictionary of observation years that most closely match each year in years
    obs_years = df_obs[year_col].unique()
    year_idx = np.arange(len(obs_years)) 
    
    #match_years = [min(years, key=lambda x:abs(x-yr)) for yr in obs_years]      
    match_years = [years[np.argmin(np.abs(np.array(years) - obs_yr))] for obs_yr in obs_years]
    year_dict = {yr: obs_years[year_idx[match_years == yr]] for yr in years}

    # Join values from the txt file that most closely matches the year of observation
    files = os.listdir(search_dir) 
    dfs = []
    columns = df_vars.index.tolist() + [index_col]
    for yr in sorted(year_dict.keys()):
        this_bn = fnmatch.filter(files, search_str % yr)[0]
        this_file = os.path.join(search_dir, this_bn)
        #import pdb; pdb.set_trace()
        df_values = pd.read_csv(this_file, sep='\t', index_col=index_col, usecols=columns)
        df_yr = df_obs[df_obs[year_col].isin(year_dict[yr])]
        # Just get records where both indices match
        df_yr = df_yr.join(df_values, how='inner')
        dfs.append(df_yr)
        #import pdb; pdb.set_trace()
        
    df_out = pd.concat(dfs)
    #import pdb; pdb.set_trace()
    print 'Length of out df: ', len(df_out)
    print 'Length of obs df: ', len(df_values)
    
    return df_out
    

def main(params):
    
    # Read params. Make variables from each line of the 1-line variables
    inputs, df_vars = read_params(params)
    '''for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])'''
    try:
        if 'years' in inputs: 
            years = np.array([int(yr) for yr in inputs['years'].split(',')])
        else:
            year_start = int(inputs['year_start'])
            year_end = int(inputs['year_end'])
            years = np.arange(year_start, year_end + 1)
        search_dir = inputs['search_dir']
        search_str = inputs['search_str']
        obs_txt = inputs['obs_txt']
        index_col = inputs['index_col']
        year_col = inputs['year_col']
        target_col = inputs['target_col']
        out_txt = inputs['out_txt']
        add_file_tag = int(inputs['add_file_tag'])
        count_type = inputs['count_type']
        day_min, day_max = [int(d) for d in inputs['day_minmax'].split(',')]
    except KeyError as e:
        missing_var = str(e).split("'")[1]
        if missing_var in ['year_start', 'year_end', 'years']:
            msg = ('No list of years or year_start/year_end specified in' +\
            ' param file:\n%s\n. Re-run script with either of these' +\
            ' parameters given.') % params
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    out_dir, out_bn = os.path.split(out_txt)
    # Add informative tags to output dir and basename
    if add_file_tag:
        res = years[1] - years[0]
        out_dir, out_dirname = os.path.split(out_dir)
        out_dirname += '_res%syr' % res
        out_dir = os.path.join(out_dir, out_dirname)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_bn = '{0}_{1}'.format(os.path.basename(obs_txt).replace('.txt', ''), out_bn)
        out_txt = os.path.join(out_dir, out_bn)
    if params != os.path.exists(os.path.join(out_dir, os.path.basename(params))):
        print 'Copying params to output dir: %s\n' % out_dir
        shutil.copy2(params, out_dir)
    
    df = get_predictors(years, search_dir, search_str, obs_txt, index_col,
                        year_col, df_vars)
    
    # Select count type and date range
    df = df[df.COUNT_TYPE == count_type]
    df.drop(['COUNT_TYPE'], axis=1, inplace=True)
    df = df[(df.DAY >= day_min) & (df.DAY <= day_max)]
    df[target_col] *= 100 # To be able to keep stuff as 8 bit ints
    df.to_csv(out_txt, sep='\t')
    
    print '\nText file written to: ', out_txt


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))
        
    