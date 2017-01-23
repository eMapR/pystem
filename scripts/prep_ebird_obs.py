# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 21:35:22 2016

@author: shooper
"""

import os
import sys
import shutil
import fnmatch
import random
import pandas as pd
import numpy as np
from osgeo import gdal
from scipy import stats
import time
import datetime
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
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))
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
    
    

def get_predictors(years, search_dir, search_str, df_obs, index_col, year_col, df_vars):
    
    # Get a dictionary of observation years that most closely match each year in years
    obs_years = df_obs[year_col].unique()
    year_idx = np.arange(len(obs_years)) 
    
    #match_years = [min(years, key=lambda x:abs(x-yr)) for yr in obs_years]      
    match_years = [years[np.argmin(np.abs(np.array(years) - obs_yr))] for obs_yr in obs_years]
    year_dict = {yr: obs_years[year_idx[match_years == yr]] for yr in years}

    # Join values from the txt file that most closely matches the year of observation
    files = os.listdir(search_dir) 
    dfs = []
    extract_cols = [c for c in df_vars.index.tolist()] + [index_col]
    for yr in sorted(year_dict.keys()):
        this_bn = fnmatch.filter(files, search_str % yr)[0]
        this_file = os.path.join(search_dir, this_bn)
        #import pdb; pdb.set_trace()
        df_values = pd.read_csv(this_file, sep='\t', index_col=index_col, usecols=extract_cols)
        df_yr = df_obs[df_obs[year_col].isin(year_dict[yr])]        
        # Just get records where both indices match
        df_yr = df_yr.join(df_values, how='inner')
        dfs.append(df_yr)
        
    df_out = pd.concat(dfs)
    
    return df_out
    

def gaussain_weights(df, val_col, max_dist, height=1):
    
    '''   f(x) = ae ** -((x - b)**2/(2c ** 2)), where a adjusts the height of the 
    curve, b adjusts position along the x axis, and c adjusts the width (stdv)
    of the curve.'''
    
    # Define gaussian function for weights
    #  x is distance from pt, a is val of pt, b is always 0 so not in equ., and c
    #  is sqrt(max_dist)
    func = lambda x, a, c: (a * np.e) ** -((np.float_(x))**2/(2 * c ** 2))
    sigma = max_dist ** .5
    for i, pt in df[[val_col,'row','col']].iterrows():
        val, row, col = pt
        # Find surrounding points
        min_row = row - max_dist/2
        max_row = row + max_dist/2
        min_col = col - max_dist/2
        max_col = col + max_dist/2
        neighbors = df[(df.row >= min_row) &
                        (df.row <= max_row) &
                        (df.col >= min_col) &
                        (df.col <= max_col)]
        if len(neighbors) == 0:
            df.ix[i, 'weighted'] = val
            continue
            
        # Calc distances from the point
        distances = np.sqrt((row - neighbors.row) ** 2 + (col - neighbors.col) ** 2)
        
        # Calc gaussian weights
        diff = val - neighbors[val_col]
        weights = func(distances, height, sigma)
        weighted_diff = np.mean(diff * weights)
        df.ix[i, 'weighted'] = val - weighted_diff
        #import pdb; pdb.set_trace()
    
    return df


def main(params, pct_train=None, aggregate_presence=False):
    t0 = time.time()
    
    # Read params. Make variables from each line of the 1-line variables
    inputs, df_vars = stem.read_params(params)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    try:
        if 'years' in inputs: 
            years = np.array([int(yr) for yr in years.split(',')])
        else:
            year_start = int(year_start)
            year_end = int(year_end)
            years = np.arange(year_start, year_end + 1)
        '''tsa_mosaic = inputs['tsa_mosaic']
        search_dir = inputs['search_dir']
        search_str = inputs['search_str']
        obs_txt = inputs['obs_txt']
        index_col = inputs['index_col']
        year_col = inputs['year_col']
        target_col = inputs['target_col']
        out_txt = inputs['out_txt']'''
        add_file_tag = int(add_file_tag)
        #count_type = inputs['count_type']
        
    except KeyError as e:
        missing_var = str(e).split("'")[1]
        if missing_var in ['year_start', 'year_end', 'years']:
            msg = ('No list of years or year_start/year_end specified in' +\
            ' param file:\n%s\n. Re-run script with either of these' +\
            ' parameters given.') % params
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    out_dir, original_bn = os.path.split(out_txt)
    # Add informative tags to output dir and basename
    if add_file_tag:
        res = years[1] - years[0]
        #out_dir = os.path.basename(out_dir)
        now = datetime.datetime.now()
        date_str = str(now.date()).replace('-','')
        time_str = str(now.time()).replace(':','')[:4]
        out_dirname = '{0}_res{1}yr_{2}_{3}'.format(target_col, res, date_str, time_str)
        out_dir = os.path.join(out_dir, out_dirname)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_bn = '{0}_{1}'.format(os.path.basename(obs_txt).replace('.txt', ''), original_bn)
        out_txt = os.path.join(out_dir, out_bn)
        
    if params != os.path.exists(os.path.join(out_dir, os.path.basename(params))):
        print 'Copying params to output dir: %s\n' % out_dir
        shutil.copy2(params, out_dir)
    
    print 'Getting predictors... '
    t1 = time.time()
    df_obs = pd.read_csv(obs_txt, sep='\t', index_col=index_col)
    original_columns = df_obs.columns
    df = get_predictors(years, search_dir, search_str, df_obs, index_col, 
                        year_col, df_vars)
    print '%.1f seconds\n' % (time.time() - t1)
    
    # Select count type and date range
    if 'count_type' in inputs:
        count_type = [t.strip() for t in count_type.split(',')]
        df = df[df.COUNT_TYPE.isin(count_type)]
        #df.drop(['COUNT_TYPE'], axis=1, inplace=True)
        if 'P21' in count_type:
            df = df[df.EFFORT_DISTANCE_KM < .1]
    if 'day_minmax' in inputs:
        day_min,  day_max  = [int(d) for d in day_minmax.split(',')]
        df = df[(df.DAY >= day_min) & (df.DAY <= day_max)]
    if 'time_minmax' in inputs:
        time_min, time_max = [int(t) for t in time_minmax.split(',')]
        df = df[(df.TIME >= time_min) & (df.TIME <= time_max)]
    if 'max_effort_time' in inputs:
        max_effort_time = int(max_effort_time)
        df = df[df.EFFORT_HRS < max_effort_time]
    if 'max_effort_dist' in inputs:
        max_effort_dist = int(max_effort_dist)
        df = df[df.EFFORT_DISTANCE_KM < max_effort_time]
    
    #df = df[(df.YEAR >= min(years)) & (df.YEAR <= max(years))]
    #df[target_col] *= 100 # To be able to keep stuff as 8 bit ints
    
    # Calc row and col from x, y
    ds = gdal.Open(tsa_mosaic)
    tx = ds.GetGeoTransform()
    ds = None
    ul_xy = tx[0], tx[3]
    df['row'], df['col']= zip(*[stem.calc_offset(ul_xy, xy, tx)
                        for i, xy in df[['x','y']].iterrows()])
    
    if 'kernel_dist' in inputs:
        t1 = time.time()
        print 'Calculating kernel density...'
        kernel_dist = int(kernel_dist)
        for yr in df.YEAR.unique():
            yr_mask = df.YEAR == yr
            df_w = gaussain_weights(df[yr_mask], target_col, kernel_dist)
            df.ix[yr_mask, target_col] = df_w.weighted
        '''
        df_w = gaussain_weights(df, target_col, kernel_dist)
        df[target_col] = df_w.weighted
        #df = df.drop_duplicates(subset=[target_col, 'row', 'col'])'''
        print '%.1f seconds\n' % (time.time() - t1)#"""
    
    if aggregate_presence:
        t1 = time.time()
        print 'Aggregating presence records...'
        df.ix[df[target_col] > 0, target_col] = 1
        for yr in df.YEAR.unique():
            yr_mask = df.YEAR == yr
            df_yr = df[yr_mask]
            # Get unique locations for this year
            unique = df_yr[['row','col']].drop_duplicates().values
            for row, col in unique:
                this_loc = df_yr[(df_yr.row == row) & (df_yr.col == col)]
                #If there are ones and 0s, drop the 0s
                if this_loc[target_col].min() == 0 and this_loc[target_col].max() == 1:
                    df.drop(this_loc[this_loc[target_col] == 0].index, inplace=True)            
        print '%.1f seconds\n' % (time.time() - t1)
        
    if pct_train:
        print 'Splitting training and test sets...'
        pct_train = float(pct_train)
        #n_test = int(len(df) * (1 - pct_train))
        unique = df[['row','col']].drop_duplicates().values
        n_test = int(len(unique) * (1 - pct_train))
        random_idx = random.sample(xrange(len(unique)), n_test)
        random_row, random_col = zip(*unique[random_idx])
        df_test = df[df.row.isin(random_row) & df.col.isin(random_col)]
        test_idx = df_test.index
        test_txt = out_txt.replace('.txt', '_test.txt')
        df = df[~df.index.isin(test_idx)]
        df_test.to_csv(test_txt, sep='\t')
    
    df.to_csv(out_txt, sep='\t')
    obs_out_txt = out_txt.replace('_' + original_bn[:-4], '')
    df[original_columns].to_csv(obs_out_txt, sep='\t')
    
    
    print '\nLength of output df:', len(df)
    print 'Text file written to: ', out_txt
    print '\nTotal time: %.1f minutes' % ((time.time() - t0)/60)


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))


"""
    if aggregate_presence:
        t1 = time.time()
        print 'Aggregating presence records...'
        '''presence_locs = df.ix[df[target_col] == 1, ['YEAR', 'row','col']].drop_duplicates()
        df_unique = df.drop_duplicates(['YEAR', 'row', 'col'])
        for i, r in presence_locs.iterrows():
            df_unique.ix[(df_unique.YEAR == r.YEAR) &
                            (df_unique.row == r.row) &
                            (df_unique.col == r.col), target_col] = 1
        df = df_unique'''
        #dfs = []
        df.ix[df[target_col] > 0, target_col] = 1
        for yr in df.YEAR.unique():
            yr_mask = df.YEAR == yr
            df_yr = df[yr_mask]
            # Get unique locations for this year
            unique = df_yr[['row','col']].drop_duplicates().values
            for row, col in unique:
                this_loc = df_yr[(df_yr.row == row) & (df_yr.col == col)]
                #If there are 1s drop all but the first (first is arbitrary)
                presence_idx = this_loc[this_loc[target_col] == 1].index
                if len(presence_idx) > 1:
                    first = presence_idx.min()
                    df.drop(presence_idx[presence_idx != first], inplace=True)
                # If there are also 0s, drop all of them
                if len(presence_idx) > 0 and this_loc[target_col].min() == 0:
                    df.drop(this_loc[this_loc[target_col] == 0].index, inplace=True)
                    import pdb; pdb.set_trace()
"""
    