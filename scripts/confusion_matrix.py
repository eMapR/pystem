# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 19:08:09 2016

@author: shooper
"""

import gdal
import os
import sys
import pandas as pd

from get_stratified_random_pixels import parse_bins
from evaluate_cover import confusion_matrix_by_area


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

    # Set the dictionary key to whatever is left of the ";" and the val
    #   to whatever is to the right. Strip whitespace too.
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))

    print '\nParameters read from:\n', txt, '\n'
    
    return d
    

def main(params, ar_p=None, out_txt=None, inventory_txt=None, target_col=None, match=False):
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
        str_check = t_path, sample_txt#, target_col
    except NameError as e:
        print ''
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
        
    # If p_path was specified, this call of the function is coming from outside
    #   predict_stem.py. Otherwise, ar_p should be given.
    if 'p_path' in locals():
        print 'Reading in the prediction path:%s\n' % p_path
        ds_p = gdal.Open(p_path)
        ar_p = ds_p.ReadAsArray()
    ds_t = gdal.Open(t_path)
    print t_path
    ar_t = ds_t.ReadAsArray()

    mask = (ar_p == p_nodata) | (ar_t == t_nodata)#'''
    
    samples = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    
    #if out_dir_: # then out_dir came from predict_stem call
    #    out_dir = out_dir_
    #out_txt = os.path.join(out_dir, 'confusion.txt')
    if out_txt:
        out_dir = os.path.dirname(out_txt)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
    
    df = confusion_matrix_by_area(ar_p, ar_t, samples, p_nodata, t_nodata, mask=mask, bins=bins, out_txt=out_txt, target_col=target_col, match=match)

    ar_p = None
    ar_t = None
    mask = None
    
    if inventory_txt:
        df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
        if 'vote' in os.path.basename(out_dir):
            df_inv.ix[s]
            ''' get stamp and save accuracy in inventory'''
    
    return df


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))