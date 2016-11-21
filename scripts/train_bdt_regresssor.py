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

def main(params, min_oob=0, err_threshold=10):
    t0 = time.time()
    
    #read_params(params)
    inputs, df_var = stem.read_params(params)

    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    try:
        if 'err_threshold' in inputs: err_threshold = int(err_threshold)
        str_check = sample_txt, target_col, mosaic_path, tsa_txt, out_dir
        n_sets = int(n_sets)
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    '''now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    stamp = '{0}_{1}_{2}'.format(target_col, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it already exists'''
    #out_dir = '/vol/v2/stem/imperv/imperv_bdt'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    stamp = os.path.basename(out_dir)
    shutil.copy2(params, out_dir) #Copy the params for reference
    
    df_train = pd.read_csv(sample_txt, sep='\t')
    # Check that df_train has exactly the same columns as variables specified in df_vars
    unmatched_vars = [v for v in df_var.index if v not in [c for c  in df_train]]
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in sample_txt but specified in params:\n' + unmatched_str
        raise NameError(msg)
    predict_cols = sorted(np.unique([c for c in df_train.columns for v in df_var.index if v in c]))
    df_var = df_var.reindex(df_var.index.sort_values())# Make sure predict_cols and df_var are in the same order
    
    # Make dataframe of set coords
    min_x, min_y, max_x, max_y, x_res, y_res, tx = stem.get_raster_bounds(mosaic_path)
    if x_res < 0:
        ul_x = max_x
        lr_x = min_x
    else:
        ul_x = min_x
        lr_x = max_x
    if y_res < 0:
        ul_y = max_y
        lr_y = min_y
    else:
        ul_y = min_y
        lr_y = max_y
    ar_sets = np.tile([ul_x, ul_y, lr_x, lr_y], n_sets).reshape(n_sets, 4)
    df_sets = pd.DataFrame(ar_sets, columns=['ul_x', 'ul_y', 'lr_x', 'lr_y'])

    # Train a tree for each support set
    print 'Training models...'
    t1 = time.time()
    #set_txt = os.path.join(out_dir, 'decisiontree_models/%s_support_sets.txt' % stamp)
    #df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    x_train = df_train.reindex(columns=predict_cols)
    y_train = df_train[target_col]
    df_sets['dt_model'] = ''
    df_sets['oob_rate'] = 0
    df_sets[['dt_model', 'oob_rate']] = [stem.fit_bdt_tree_regressor(x_train, y_train) for s in df_sets.index]
    del df_train
    print 'Estimated average OOB score: ', int(df_sets.oob_rate.mean())
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    # Write df_sets and each decison tree to disk
    print 'Saving models...'
    t1 = time.time()
    df_sets, set_txt = stem.write_model(out_dir, df_sets)
    print '%.1f minutes\n' % ((time.time() - t1)/60)#'''
    
    '''out_dir = '/vol/v2/stem/canopy/models/canopy_20161016_0910'
    stamp = os.path.basename(out_dir)
    set_txt = '/vol/v2/stem/{0}/models/{1}/decisiontree_models/{1}_support_sets.txt'.format(target_col, stamp)
    df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    oob_txt = os.path.join(out_dir, '%s_oob.txt' % stamp)
    df_oob = pd.read_csv(oob_txt, sep='\t')
    ds = gdal.Open(os.path.join(out_dir, '%s_oob.bsq' % stamp))
    ar_oob = ds.ReadAsArray()
    ds = None
    ds = gdal.Open(os.path.join(out_dir, '%s_count.bsq' % stamp))
    ar_cnt = ds.ReadAsArray()
    ds = None 
    
    predict_cols = ['aspectNESW','aspectNWSE','brightness','delta_bright','delta_green','delta_nbr','delta_wet', 'elevation','greenness','mse','nbr','slope','time_since','wetness']#'''
    
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
        n_samples = int(sample_txt.split('_')[1].replace('sample',''))
        info = [model_type, None, None, None, None, None, None, None, None, n_sets, n_samples, str(support_size), sets_per_cell, max_features]
        df_inv.ix[stamp] = info
        info_dir = os.path.dirname(inventory_txt)
        existing_models = fnmatch.filter(os.listdir(os.path.dirname(info_dir)), '%s*' % target_col)
        if len(existing_models) > 0:
            df_inv = df_inv[df_inv.index.isin(existing_models)]
        
        # Check if oob_map params were specified. If not, set to defaults
        if 'err_threshold' not in locals():
            print 'err_threshold not specified. Using default: 10 ...\n'
            err_threshold = 10
        else:
            err_threshold = int(err_threshold)
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
                                     predict_cols, err_threshold, out_dir,
                                     stamp, prj, driver)
        df_sets.to_csv(set_txt, sep='\t')#'''
    
        #if 'inventory_txt' in locals() :
        avg_oob = round(np.mean(ar_oob[mask]), 1)
        avg_cnt = int(round(np.mean(ar_cnt[mask]), 0))
        df_inv.ix[stamp, 'avg_oob'] = avg_oob
        df_inv.ix[stamp, 'avg_count'] = avg_cnt
        if len(df_inv) > 1:
            df_inv.to_csv(inventory_txt, sep='\t')
        else:
            print 'WARNING: Model info not written to inventory_txt...\n'
        
        print '\nAverage OOB score: .................... %.1f' % avg_oob
        print '\nAverage number of overlapping sets: ... %s\n' % avg_cnt
        
        print 'Time to make OOB score map: %.1f hours\n' % ((time.time() - t1)/3600)
            
        #except Exception as e:
        #    print 'Problem getting oob map: ', e
    
    print 'Total training time: %.1f minutes' % ((time.time() - t0)/60)

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))
