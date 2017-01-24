# -*- coding: utf-8 -*-
"""
Generate stratified samples of specified bin sizes from a given input raster 

@author: Sam Hooper, samhooperstudio@gmail.com
"""
import gdal
import random
import sys
import os
import time
import shutil
import pandas as pd
import numpy as np
from datetime import datetime

import extract_xy_by_tsa as extract

def read_params(txt):
    '''
    Return a dictionary from parsed parameters in txt
    '''
    d = {}
    
    # Read in the text file
    try:
        with open(txt) as f:
            input_vars = [line.split(";") for line in f]
    except: 
        print 'Problem reading parameter file: ', txt
        return None
    
    # Set the dictionary key to whatever is left of the ";" and the val
    #   to whatever is to the right. Strip whitespace too.
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))
        
    print 'Parameters read from:\n', txt, '\n'
    return d


def parse_bins(bin_str):
    ''' Integerize a string of min:max and return as a list of length 2'''
    bin_list = [b.split(':') for b in bin_str.split(',')]
    bins = [(int(mn), int(mx)) for mn, mx in bin_list]
    
    return bins


def extract_by_kernel(ar, rows, cols, data_type, col_name, nodata):
    
    row_dirs = [-1,-1,-1, 0, 0, 0, 1, 1, 1]
    col_dirs = [-1, 0, 1,-1, 0, 1,-1, 0, 1]
    kernel_rows = [row + d + 1 for row in rows for d in row_dirs]
    kernel_cols = [col + d + 1 for col in cols for d in col_dirs]
    
    ar_buf = np.full([dim + 2 for dim in ar.shape], nodata, dtype=np.int32)
    ar_buf[1:-1, 1:-1] = ar
    del ar
    
    kernel_vals = ar_buf[kernel_rows, kernel_cols].reshape(len(rows), len(row_dirs))
    train_stats = pd.DataFrame(extract.calc_row_stats(kernel_vals, data_type, col_name, nodata))
    vals = train_stats[col_name].astype(np.int32)
    
    return vals
    

def get_stratified_sample(raster_path, col_name, data_band, n_samples, bins, pct_train=None, nodata=None, zero_inflation=None, data_type='continuous', kernel=False):
    '''
    Return a dataframe of stratified randomly sampled pixels from raster_path
    '''
    print 'Reading the raster_path... %s\n' % datetime.now()
    ds = gdal.Open(raster_path)
    tx = ds.GetGeoTransform()
    band = ds.GetRasterBand(data_band)
    #ar_data = band.ReadAsArray()
    ar = band.ReadAsArray()
    if nodata == None:
        nodata = band.GetNoDataValue()
        if nodata == None:
            sys.exit('Could not obtain nodata value from dataset and' +\
            ' none specified in parameters file. Try re-running with' +\
            'nodata specified.')
    
    # Buffer the array by 1 on each side to be able to extract 3x3 kernels
    #ar = np.full([dim + 2 for dim in ar_data.shape], nodata, dtype=np.int32)
    #ar[1:-1, 1:-1] = ar_data
    
    samples_per = n_samples/len(bins)
    
    print 'Making arrays of row and col indices... %s\n' % datetime.now()
    # Make 2D arrays of row and column indices
    ar_rows, ar_cols = np.indices(ar.shape)
    #ar_rows = ar_rows[1:-1, 1:-1] # Only sample inds that have data
    #ar_cols = ar_cols[1:-1, 1:-1]

    
    # Get random samples for each bin
    train_rows = []
    train_cols = []
    test_rows = []
    test_cols = []
    #nodata_mask = np.logical_and(ar != nodata, ind_mask)
    #nodata_mask = ar_data != nodata
    nodata_mask = ar != nodata
    for b in bins:
        t1 = time.time()
        this_min, this_max = b
        print 'Getting random samples between %s and %s...' % (this_min, this_max)
        mask = (ar > this_min) & (ar <= this_max) & nodata_mask
        these_rows = ar_rows[mask]
        these_cols = ar_cols[mask]
        
        try:
            if this_max == 0 and zero_inflation:
                samples = random.sample(xrange(len(these_rows)), samples_per * zero_inflation)
            else:
                samples = random.sample(xrange(len(these_rows)), samples_per)
            tr_rows = these_rows[samples]
            tr_cols = these_cols[samples]
        # If there aren't enough pixels to generate samples for this bin
        except:
            print ('Not enough pixels between %s and %s to generate %s' +\
            ' random samples. Returning all pixels for this bin.')\
            % (this_min, this_max, samples_per)
            tr_rows = these_rows
            tr_cols = these_cols
        
        # If pct_train is specified, split the sample indices into train/test sets
        te_rows = []
        te_cols = []
        if pct_train:
            split_ind = len(tr_rows) * pct_train
            te_rows = tr_rows[split_ind:]
            te_cols = tr_cols[split_ind:]
            tr_rows = tr_rows[:split_ind]
            tr_cols = tr_cols[:split_ind]

        train_rows.extend(tr_rows)
        train_cols.extend(tr_cols)
        test_rows.extend(te_rows)
        test_cols.extend(te_cols)
        print '%.1f seconds\n' % (time.time() - t1)
    del tr_rows, tr_cols, te_rows, te_cols
    
    # If True, extract with 3x3 kernel. Otherwise, just get the vals (row,col)
    if kernel:
        train_vals = extract_by_kernel(ar, train_rows, train_cols, data_type, col_name, nodata)
    else:
        train_vals = ar[train_rows, train_cols]
        
    # Calculate x and y for later extractions
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    train_x = [int(ul_x + c * x_res) for c in train_cols]
    train_y = [int(ul_y + r * y_res) for r in train_rows]
    df_train = pd.DataFrame({'x': train_x,
                             'y': train_y,
                             'row': train_rows,
                             'col': train_cols,
                             col_name: train_vals
                             })
                             
    # If training and testing samples were split, get test vals                    
    df_test = None
    if pct_train:
        if kernel:
            test_vals = extract_by_kernel(ar, train_rows, train_cols, data_type, col_name, nodata)
        else:
            test_vals = ar[test_rows, test_cols]
        test_x = [int(ul_x + c * x_res) for c in test_cols]
        test_y = [int(ul_y + r * y_res) for r in test_rows]
        df_test = pd.DataFrame({'x': test_x, 
                                'y': test_y, 
                                'row': test_rows,
                                'col': test_cols,
                                col_name: test_vals
                                })
    
    return df_train, df_test


def main(params, data_band=1, nodata=None, data_type='continuous', kernel=False):
    
    t0 = time.time()
    data_band = None
    nodata = None
    zero_inflation = None
    
    # Read params and make variables from each line
    inputs = read_params(params)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    
    out_dir = os.path.dirname(out_txt)
    '''if not os.path.exists(out_dir):
        print 'Warning: output directory does not exist. Creating directory...'
        os.makedirs(out_dir)'''
    
    # Integerize numeric params
    if 'data_band' in locals(): data_band = int(data_band)
    if 'nodata' in locals(): nodata = int(nodata)
    if 'pct_train' in locals(): 
        pct_train = float(pct_train)
    else: 
        pct_train = None
    if zero_inflation: zero_inflation = int(zero_inflation)
    n_samples = int(n_samples)
    bins = parse_bins(bins)
    
    # Generate samples
    df_train, df_test = get_stratified_sample(raster_path, col_name, data_band,
                                              n_samples, bins, pct_train, 
                                              nodata, zero_inflation, data_type,
                                              kernel)
    df_train['obs_id'] = df_train.index
    
    # Write samples to text file
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    stamp = '{0}_{1}_{2}'.format(len(df_train), date_str, time_str)
    out_txt = out_txt.replace('.txt', '%s.txt' % stamp)
    bn = os.path.basename(out_txt)
    out_dir = os.path.join(os.path.dirname(out_txt), bn[:-4])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_txt = os.path.join(out_dir, bn)
    df_train.to_csv(out_txt, sep='\t', index=False)
    print 'Samples written to:\n%s\n' % out_txt
    
    shutil.copy2(params, out_dir) #Copy the params for reference
    
    if pct_train:
        df_test['obs_id'] = df_test.index
        test_txt = out_txt.replace('%s.txt' % stamp, '%s_test.txt' % stamp)
        df_test.to_csv(test_txt, sep='\t', index=False)
        print 'Test samples written to directory:\n%s' % out_dir
    
    print 'Total time: %.1f minutes' % ((time.time() - t0)/60)
        
        
if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))

#raster_path = '/vol/v2/stem/canopy/truth_map/nlcd_canopy_2001_nodata.tif'
#out_txt = '/home/server/student/homes/shooper/scripts/random_pixel_test.txt'
#params = '/vol/v2/stem/scripts/get_stratified_random_pixels_params.txt'
#df = main(params)