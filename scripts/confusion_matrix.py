# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 19:08:09 2016

@author: shooper
"""

import gdal
import os
import sys
import shutil
import pandas as pd
import numpy as np

from get_stratified_random_pixels import parse_bins
from evaluation import confusion_matrix_by_area
import mosaic_by_tsa as mosaic


def read_params(txt):
    '''
    Return a dictionary and a dataframe from parsed parameters in txt
    '''
    if not os.path.exists(txt):
        print 'Param file does not exist:\n%s' % txt
        return None

    # Read in the rest of the text file line by line
    d = {}
    try:
        with open(txt) as f:
            input_vars = [line.split(";") for line in f]
    except:
        print 'Problem reading parameter file:\n', txt
        return None

    # Set the dictionary key to whatever is left of the ";" and the value
    #   to whatever is to the right. Strip whitespace too.
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))

    print '\nParameters read from:\n', txt, '\n'
    
    return d
    

def main(params, ar_p=None, out_txt=None, inventory_txt=None, target_col=None, match=False, file_stamp=None):
    #p_path, t_path, bins, sample_txt, p_nodata, t_nodata, out_dir, inventory_txt=None
    
    # Read params and make variables from text
    inputs = read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])

    # Check that variables were specified in params
    try:
        bins = parse_bins(bins)
        p_nodata = int(p_nodata)
        t_nodata = int(t_nodata)
        str_check = sample_txt#, target_col
    except NameError as e:
        print ''
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    #if out_dir_: # then out_dir came from predict_stem call
    #    out_dir = out_dir_
    #out_txt = os.path.join(out_dir, 'confusion.txt')
    if out_txt:
        out_dir = os.path.dirname(out_txt)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)    
        shutil.copy2(params, out_dir)
        
    # If p_path was specified, this call of the function is coming from outside
    #   predict_stem.py. Otherwise, ar_p should be given.
    if 'p_path' in locals():
        print 'Reading in the prediction raster:%s\n' % p_path
        ds_p = gdal.Open(p_path)
        ar_p = ds_p.ReadAsArray()
        
    
    ds_t = gdal.Open(t_path)
    ar_t = ds_t.ReadAsArray()
    
    t_xsize = ds_t.RasterXSize
    t_ysize = ds_t.RasterYSize
    p_xsize = ds_p.RasterXSize
    p_ysize = ds_p.RasterYSize
    
    # If two arrays are different sizes, make prediction array match reference
    if not t_xsize == p_xsize or t_ysize == p_ysize:
        tx_t = ds_t.GetGeoTransform()
        tx_p = ds_p.GetGeoTransform()
        offset = mosaic.calc_offset((tx_t[0], tx_t[3]), tx_p)
        t_inds, p_inds = mosaic.get_offset_array_indices((t_ysize, t_xsize), (p_ysize, p_xsize), offset)
        ar_buf = np.full(ar_t.shape, p_nodata, dtype=ar_p.dtype)
        ar_buf[t_inds[0]:t_inds[1], t_inds[2]:t_inds[3]] = ar_p[p_inds[0]:p_inds[1], p_inds[2]:p_inds[3]]
        ar_p = ar_buf.copy()
        del ar_buf
    mask = (ar_p == p_nodata) | (ar_t == t_nodata)#'''
    
    samples = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    
    df = confusion_matrix_by_area(ar_p, ar_t, samples, p_nodata, t_nodata, mask=mask, bins=bins, out_txt=out_txt, target_col=target_col, match=match)

    ar_p = None
    ar_t = None
    mask = None
    
    accuracy = df.ix['producer','user']
    kappa = df.ix['producer', 'kappa']
    if inventory_txt and file_stamp:
        df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
        if file_stamp in df_inv.index and 'vote' in os.path.basename(out_dir):
            cols = ['vote_accuracy', 'vote_kappa']
            df_inv.ix[file_stamp, cols] = accuracy, kappa
            df_inv.to_csv(inventory_txt, sep='\t')
            print 'Vote scores written to inventory_txt: ', inventory_txt
            
        if file_stamp in df_inv.index and 'mean' in os.path.basename(out_dir):
            cols = ['mean_accuracy', 'mean_kappa']
            df_inv.ix[file_stamp, cols] = accuracy, kappa
            df_inv.to_csv(inventory_txt, sep='\t')
            
    
    return df


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))