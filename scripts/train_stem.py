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
import pandas as pd
import numpy as np
from datetime import datetime

# Import ancillary scripts
import stem
import generate_gsrd as gsrd

def main(params):
    
    t0 = time.time()
    #read_params(params)
    inputs, df_var = stem.read_params(params)

    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    try:
        if 'max_features' not in locals(): max_features=None
        if 'pct_train' not in locals(): pct_train=None
        num_vars = stem.vars_to_numbers(cell_size, support_size, sets_per_cell,
                                   min_obs, max_features, pct_train)
        cell_size, support_size, sets_per_cell, min_obs, max_features, pct_train = num_vars
        str_check = sample_txt, target_col, mosaic_path, tsa_txt, out_dir
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    stamp = '{0}_{1}_{2}'.format(target_col, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it already exists
    shutil.copy2(params, out_dir) #Copy the params for reference
    
    # Get samples and support set bounds
    if 'gsrd_shp' not in locals(): gsrd_shp = None
    out_txt = os.path.join(out_dir, stamp + '.txt')
    df_train, df_sets = gsrd.get_gsrd(mosaic_path, cell_size, support_size, sets_per_cell,
                        sample_txt, min_obs,  target_col, out_txt,
                        gsrd_shp, pct_train)

    # Check that df_train has exactly the same columns as variables specified in df_vars
    unmatched_vars = [v for v in df_var.index if v not in [c for c  in df_train]]
    
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in sample_txt but specified in params:\n' + unmatched_str
        raise NameError(msg)
    
    predict_cols = sorted(np.unique([c for c in df_train.columns for v in df_var.index if v in c]))
    df_var = df_var.reindex(df_var.index.sort_values())# Make sure predict_cols and df_var are in the same order

    # Train a tree for each support set
    x_train = df_train.reindex(columns=predict_cols + ['set_id'])
    y_train = df_train[[target_col, 'set_id']]    
    df_sets['dt_model'] = [stem.fit_tree(x_train.ix[x_train.set_id==s, predict_cols],\
    y_train.ix[y_train.set_id==s, target_col], max_features) for s in df_sets.index]
    
    # Write df_sets and each decison tree to disk
    stem.write_model(out_dir, df_sets)
    
    # Record params in inventory text file
    if 'inventory_txt' in locals():
        t1 = time.time()
        print 'Getting model info... '
        df_inv = pd.read_csv(inventory_txt, sep='\t')
        if 'regressor' in params: 
            model_type = 'Regressor'
        else: 
            model_type = 'Classifier'
        n_sets = len(df_sets)
        n_samples = int(sample_txt.split('_')[1].replace('sample',''))
        info = [stamp, model_type, '', '', '', '', '','', n_sets, n_samples, str(support_size), sets_per_cell, max_features]
        df_inv = df_inv.append(pd.DataFrame([info], columns=df_inv.columns),ignore_index=True)
        existing_models = fnmatch.filter(os.listdir(os.path.dirname(out_dir)), '%s*' % target_col)
        df_inv = df_inv[df_inv.stamp.isin(existing_models)]
        df_inv.to_csv(inventory_txt, sep='\t', index=False)
        print '%s seconds\n' % (time.time() - t1)
    
    if 'oob_txt' in locals() and 'mosaic_path' in locals():        
        t1 = time.time()
        print 'Calculating OOB scores and making OOB map...'
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
        
        df_oob = pd.DataFrame(oob_txt, sep='\t')
        import get_oob_map as oob
        ar_oob, ar_cnt = oob.oob_map(ysize, xsize, nodata, mask, n_tiles, tx, support_size, df_oob, df_sets, val_col, var_cols, err_threshold, out_dir, file_stamp, prj, driver)
    
        if 'inventory_txt' in locals() :
            avg_oob = np.mean(ar_oob[~mask])
            avg_cnt = int(round(np.mean(ar_cnt[~mask]), 0))
            df_inv.ix[file_stamp, 'oob_score'] = avg_oob
            df_inv.ix[file_stamp, 'avg_nsets'] = avg_cnt
            print '\nAverage OOB score: .................... %.1f' % avg_oob
            print '\nAverage number of overlapping sets: ... %s\n' % avg_count
        
        print '%.1f minutes' % ((time.time() - t1)/60)
        
    
    print 'Total training time: %.1f minutes' % ((time.time() - t0)/60)

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))
