# -*- coding: utf-8 -*-
"""
Train and write to disk a spatiotemporal exploratory model
 

@author: shooper
"""
import os
import sys
import time
import gdal
import shutil
import fnmatch
from datetime import datetime
import pandas as pd
import numpy as np

# Import ancillary scripts
import stem
#import generate_gsrd as gsrd

def main(params, pct_train=None, min_oob=0):
    t0 = time.time()
    
    #read_params(params)
    inputs, df_var = stem.read_params(params)
    
    # Convert params to named variables and check for required vars
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    try:
        if 'max_features' not in locals(): max_features=None
        if 'min_oob' in inputs: min_oob = int(min_oob)
        num_vars = stem.vars_to_numbers(cell_size, support_size, sets_per_cell,
                                   min_obs, max_features, pct_train)
        cell_size, support_size, sets_per_cell, min_obs, max_features, pct_train = num_vars
        str_check = sample_txt, target_col, mosaic_path, tsa_txt, out_dir, model_type
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    # Make a timestamped output directory if outdir not specified
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    if not 'out_dirname' in locals(): out_dirname = target_col
    stamp = '{0}_{1}_{2}'.format(out_dirname, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it already exists'''
    shutil.copy2(params, out_dir) #Copy the params for reference

    # Read in training samples and check that df_train has exactly the same
    #   columns as variables specified in df_vars    
    df_train = pd.read_csv(sample_txt, sep='\t')
    n_samples = len(df_train)
    unmatched_vars = [v for v in df_var.index if v not in [c for c  in df_train]]
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in sample_txt but specified in params:\n' + unmatched_str
        raise NameError(msg)
    predict_cols = sorted(np.unique([c for c in df_train.columns for v in df_var.index if v in c]))
    df_var = df_var.reindex(df_var.index.sort_values())# Make sure predict_cols and df_var are in the same order
    
    # If there are variables that should remain constant across the modeling
    #   region, get the names
    if 'constant_vars' in locals():
        constant_vars = sorted([i.strip() for i in constant_vars.split(',')])
        predict_cols += constant_vars
    
    # Get samples and support set bounds
    if 'gsrd_shp' not in locals(): gsrd_shp = None
    out_txt = os.path.join(out_dir, stamp + '.txt')
    df_train, df_sets, df_oob = stem.get_gsrd(mosaic_path, cell_size, support_size,
                                              sets_per_cell, df_train, min_obs,
                                              target_col, predict_cols, out_txt,
                                              gsrd_shp, pct_train)

    # Train a tree for each support set
    print 'Training models...'
    t1 = time.time()
    if model_type.lower() == 'classifier':
        print 'Training STEM with classifier algorithm...'
        model_func = stem.fit_tree_classifier
    else:
        print 'Training STEM with regressor algorithm...'
        model_func = stem.fit_tree_regressor
    x_train = df_train.reindex(columns=predict_cols + ['set_id'])
    y_train = df_train[[target_col, 'set_id']]    
    df_sets['dt_model'] = [model_func(x_train.ix[x_train.set_id==s, predict_cols],\
    y_train.ix[y_train.set_id==s, target_col], max_features) for s in df_sets.index]
    del df_train
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    # Calculate OOB rates and drop sets with too low OOB
    print 'Calculating OOB rates...'
    t1 = time.time()
    df_sets, low_oob = stem.get_oob_rates(df_sets, df_oob, target_col, predict_cols, min_oob)
    if len(low_oob) > 0:
        df_sets.drop(low_oob.index, inplace=True)
        low_oob_shp = os.path.join(out_dir, 'gsrd_low_oob.shp')
        low_oob.drop('dt_model', axis=1, inplace=True)
        stem.coords_to_shp(low_oob, gsrd_shp, low_oob_shp)
    print '%s sets dropped because OOB rate < %s' % (len(low_oob), min_oob)
    print 'Min OOB rate after dropping: ', df_sets.oob_rate.min()
    print 'Estimated average OOB score: ', int(df_sets.oob_rate.mean())
    print '%.1f minutes\n' % ((time.time() - t1)/60)

    # Write df_sets and each decison tree to disk
    print 'Saving models...'
    t1 = time.time()
    df_sets, set_txt = stem.write_model(out_dir, df_sets)
    print '%.1f minutes\n' % ((time.time() - t1)/60)#'''
    
    #stamp = os.path.basename(out_dir)
    #set_txt = '/vol/v2/stem/{0}/models/{1}/decisiontree_models/{1}_support_sets.txt'.format(target_col, stamp) 
    
    #predict_cols = ['aspectNESW','aspectNWSE','brightness','delta_bright','delta_green','delta_nbr','delta_wet', 'elevation','greenness','mse','nbr','slope','time_since','wetness']#'''
    
    # Record params in inventory text file
    if 'inventory_txt' in locals():
        t1 = time.time()
        print 'Getting model info...\n'
        df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
        if 'regressor' in params: 
            model_type = 'Regressor'
        else: 
            model_type = 'Classifier'
        n_sets = len(df_sets)
        if 'sample' in sample_txt:
            n_samples = int(sample_txt.split('_')[1].replace('sample',''))
        info = [model_type, None, None, None, None, None, None, None, None, n_sets, n_samples, str(support_size), sets_per_cell, max_features]
        df_inv.ix[stamp] = info
        info_dir = os.path.dirname(inventory_txt)
        existing_models = fnmatch.filter(os.listdir(os.path.dirname(info_dir)), '%s*' % target_col)
        if len(existing_models) > 0:
            df_inv = df_inv[df_inv.index.isin(existing_models)]#'''
        
        # Check if oob_map params were specified. If not, set to defaults
        if 'n_tiles' not in locals():
            print 'n_tiles not specified. Using default: 25 x 15 ...\n'
            n_tiles = 25, 15
        else:
            n_tiles = int(n_tiles[0]), int(n_tiles[1])
            
        #t1 = time.time()
        print 'Calculating OOB score and making OOB score map...'
        ds = gdal.Open(mosaic_path)
        ar = ds.ReadAsArray()
        mask = ar != 0
        del ar
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        tx = ds.GetGeoTransform()
        prj = ds.GetProjection()
        driver = ds.GetDriver()
        ds = None  
        
        #import get_oob_map as oob
        ar_oob, ar_cnt, df_sets = stem.oob_map(ysize, xsize, 0, mask, n_tiles, tx,
                                     support_size, df_oob, df_sets, target_col,
                                     predict_cols, out_dir,
                                     stamp, prj, driver)
        df_sets.to_csv(set_txt, sep='\t')#'''
    
        #if 'inventory_txt' in locals() :
        avg_oob = round(np.mean(ar_oob[mask]), 1)
        avg_cnt = int(round(np.mean(ar_cnt[mask]), 0))
        df_inv.ix[stamp, 'avg_oob'] = avg_oob
        #df_inv.ix[stamp, 'avg_count'] = avg_cnt
        if len(df_inv) > 1:
            df_inv.to_csv(inventory_txt, sep='\t')
        else:
            print 'WARNING: Model info not written to inventory_txt...\n' #'''
        
        print '\nAverage OOB score: .................... %.1f' % avg_oob
        print '\nAverage number of overlapping sets: ... %s\n' % avg_cnt
        
        print 'Time to make OOB score map: %.1f hours\n' % ((time.time() - t1)/3600)
        
    
    print 'Total training time: %.1f minutes' % ((time.time() - t0)/60)

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))
