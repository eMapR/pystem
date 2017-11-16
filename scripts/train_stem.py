# -*- coding: utf-8 -*-
"""
Train and write to disk a spatiotemporal exploratory model
 

@author: shooper
"""
import os
import sys
import time
import shutil
import warnings
from datetime import datetime
from osgeo import gdal, ogr
import pandas as pd
import numpy as np
import cPickle as pickle

# Import ancillary scripts
import stem
#import generate_gsrd as gsrd

def main(params, pct_train=None, min_oob=0, gsrd_shp=None, resolution=30, make_oob_map=False, snap_coord=None, oob_map_metric='oob_rate'):
    t0 = time.time()
    
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
        str_check = sample_txt, target_col, mosaic_path, out_dir, model_type
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)

    # Read in training samples and check that df_train has exactly the same
    #   columns as variables specified in df_vars    
    df_train = pd.read_csv(sample_txt, sep='\t')
    n_samples = len(df_train)
    unmatched_vars = [v for v in df_var.index if v not in [c for c  in df_train]]
    if len(unmatched_vars) != 0:
        unmatched_str = '\n\t'.join(unmatched_vars)
        msg = 'Columns not in sample_txt but specified in params:\n\t' + unmatched_str
        import pdb; pdb.set_trace()
        raise NameError(msg)
    if target_col not in df_train.columns:
        raise NameError('target_col "%s" not in sample_txt: %s' % (target_col, sample_txt))
    
    '''# Make a timestamped output directory if outdir not specified
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    if not 'out_dirname' in locals(): out_dirname = target_col
    stamp = '{0}_{1}_{2}'.format(out_dirname, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it already exists
    shutil.copy2(params, out_dir) #Copy the params for reference '''

    predict_cols = sorted(np.unique([c for c in df_train.columns for v in df_var.index if v in c]))
    df_var = df_var.reindex(df_var.index.sort_values())# Make sure predict_cols and df_var are in the same order
    
    # If there are variables that should remain constant across the modeling
    #   region, get the names
    if 'constant_vars' in locals():
        constant_vars = sorted([i.strip() for i in constant_vars.split(',')])
        predict_cols += constant_vars
    
    # Get samples and support set bounds
    if 'gsrd_shp' not in locals(): gsrd_shp = None
    if snap_coord:
        snap_coord = [int(c) for c in snap_coord.split(',')]
    out_txt = os.path.join(out_dir, stamp + '.txt')
    '''train_dict, df_sets, oob_dict = stem.get_gsrd(mosaic_path, cell_size, support_size,
                                              sets_per_cell, df_train, min_obs,
                                              target_col, predict_cols, out_txt,
                                              gsrd_shp, pct_train, snap_coord=snap_coord)#'''
    df_sets = stem.get_gsrd(mosaic_path, cell_size, support_size,
                                              sets_per_cell, df_train, min_obs,
                                              target_col, predict_cols, out_txt,
                                              gsrd_shp, pct_train, snap_coord=snap_coord)
    n_sets = len(df_sets)            
    '''out_dir = '/vol/v2/stem/conus/models/canopy_20171019_1520'
    stamp = os.path.basename(out_dir)
    oob_pkl = os.path.join(out_dir, '%s_oob_dict.pkl' % stamp)
    with open(oob_pkl, 'rb') as f:
        oob_dict = pickle.load(f)
    train_pkl = os.path.join(out_dir, '%s_train_dict.pkl' % stamp)
    with open(train_pkl, 'rb') as f:
        train_dict = pickle.load(f)
    df_sets = pd.read_csv(os.path.join(out_dir, ))#'''
    
    # Train a tree for each support set
    t1 = time.time()
    if model_type.lower() == 'classifier':
        print 'Training STEM with classifier algorithm...'
        model_func = stem.fit_tree_classifier
    else:
        print 'Training STEM with regressor algorithm...'
        model_func = stem.fit_tree_regressor
    x_train = df_train.reindex(columns=predict_cols)
    y_train = df_train[target_col]
    importance_cols = ['importance_%s' % c for c in predict_cols]
    for c in importance_cols:
        df_sets[c] = 0
    
    # Train models
    dropped_sets = pd.DataFrame(columns=df_sets.columns)
    dt_dir = os.path.join(out_dir, 'decisiontree_models')
    if not os.path.exists(dt_dir):
        os.mkdir(dt_dir)
    dt_path_template = os.path.join(dt_dir, stamp + '_decisiontree_%s.pkl')
    oob_dict = {}
    oob_rates = [0]
    for i, (ind, ss) in enumerate(df_sets.iterrows()):
        format_obj = i + 1, n_sets, float(i)/n_sets * 100, (time.time() - t1)/60, int(np.mean(oob_rates))
        sys.stdout.write('\rTraining %s of %s decision trees (%.1f%%). Cumulative time: %.1f minutes. Avg oob: %s' % format_obj)
        sys.stdout.flush()

        # Get all samples within support set
        sample_inds = df_train.index[(df_train['x'] > ss[['ul_x', 'lr_x']].min()) &
        (df_train['x'] < ss[['ul_x', 'lr_x']].max()) &
        (df_train['y'] > ss[['ul_y', 'lr_y']].min()) &
        (df_train['y'] < ss[['ul_y', 'lr_y']].max())]
        
        n_samples = int(len(sample_inds) * .63)
        if n_samples < min_obs:
            df_sets.drop(ind, inplace=True)
            continue
        
        this_x = x_train.ix[sample_inds]
        this_y = y_train.ix[sample_inds]
        support_set = df_sets.ix[ind]
        dt_path = dt_path_template % ind
        dt_model, train_inds, oob_inds, importance, oob_metrics = stem.train_estimator(support_set, n_samples, this_x, this_y, model_func, model_type, max_features, dt_path)
        oob_rates.append(oob_metrics['oob_rate'])
        df_sets.ix[ind, importance_cols] = importance
        df_sets.ix[ind, 'dt_model'] = dt_model
        df_sets.ix[ind, 'dt_file'] = dt_path
        df_sets.ix[ind, 'n_samples'] = n_samples
        for metric in oob_metrics:
            df_sets.ix[ind, metric] = oob_metrics[metric]
        
        # Save oob and train inds #### Change to database structure
        t_ind_path = dt_path.replace('.pkl', '_train_inds.txt')
        with open(t_ind_path, 'w') as f:
            [f.write('%s\n' % ti) for ti in train_inds]
        o_ind_path = dt_path.replace('.pkl', '_oob_inds.txt')
        with open(o_ind_path, 'w') as f:
            [f.write('%s\n' % oi) for oi in oob_inds]
        oob_dict[ind] = oob_inds

    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    # Calculate OOB rates and drop sets with too low OOB
    print 'Calculating OOB rates...'
    t1 = time.time()
    df_sets, low_oob = stem.get_oob_rates(df_sets, df_train, oob_dict, target_col, predict_cols, min_oob)
    if len(low_oob) > 0:
        #df_sets.drop(low_oob.index, inplace=True)
        low_oob_shp = os.path.join(out_dir, 'low_oob_sets.shp')
        low_oob.drop('dt_model', axis=1, inplace=True)
        stem.coords_to_shp(low_oob, gsrd_shp, low_oob_shp)
    set_shp = os.path.join(out_dir, 'support_sets.shp')
    try:
        stem.coords_to_shp(df_sets, gsrd_shp, set_shp)
    except Exception as e:
        print e.message
    print '%s sets dropped because OOB rate < %s' % (len(low_oob), min_oob)
    print 'Min OOB rate after dropping: ', df_sets.oob_rate.min()
    print 'Estimated average OOB score: ', int(df_sets.oob_rate.mean())
    print '%.1f minutes\n' % ((time.time() - t1)/60)

    # Write df_sets and each decison tree to disk
    print 'Saving models...'
    set_txt = os.path.join(dt_dir, stamp + '_support_sets.txt')
    df_sets['set_id'] = df_sets.index
    df_sets.drop('dt_model', axis=1).to_csv(set_txt, sep='\t', index=False)
    t1 = time.time()
    #df_sets, set_txt = stem.write_model(out_dir, df_sets)
    print '%.1f minutes\n' % ((time.time() - t1)/60) #"""
    
    stamp = os.path.basename(out_dir)
    set_txt = '/vol/v2/stem/conus/models/{1}/decisiontree_models/{1}_support_sets.txt'.format(target_col, stamp) 
    df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    oob_pkl = os.path.join(out_dir, '%s_oob_dict.pkl' % stamp)
    #with open(oob_pkl, 'rb') as f:
    #    oob_dict = pickle.load(f)
    oob_dict = {}
    predict_cols = ['aspectNESW','aspectNWSE','brightness','delta_brightness','delta_greenness','delta_nbr','delta_wetness', 'elevation','greenness','mse','nbr','slope','time_since','wetness']#'''
    if make_oob_map:
        # Check if oob_map params were specified. If not, set to defaults
        if 'n_tiles' not in locals():
            n_tiles = 90, 40
            print 'n_tiles not specified. Using default: %s x %s ...\n' % (n_tiles)
            
        else:
            n_tiles = int(n_tiles[0]), int(n_tiles[1])
            
        print 'Calculating OOB score and making OOB score map...'
        try:
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
        except:
            mosaic_ds = ogr.Open(mosaic_path)
            if 'resolution' not in inputs:
                warnings.warn('Resolution not specified. Assuming default of 30...\n')
            mask = mosaic_ds.GetLayer()
            min_x, max_x, min_y, max_y = mask.GetExtent()
            ul_x = min_x - ((min_x - snap_coord[0]) % resolution)
            ul_y = max_y - ((max_y - snap_coord[1]) % resolution)
            xsize = int((max_x - ul_x)/resolution)
            ysize = int((ul_y - min_y)/resolution)
            prj = mask.GetSpatialRef().ExportToWkt()
            driver = gdal.GetDriverByName('gtiff')
            x_res = resolution
            y_res = -resolution
            tx = ul_x, x_res, 0, ul_y, 0, y_res
        
        avg_dict, df_sets = stem.oob_map(ysize, xsize, 0, mask, n_tiles, tx,
                                     support_size, oob_dict, df_sets, df_train, target_col,
                                     predict_cols, out_dir,
                                     stamp, prj, driver, oob_map_metric)
        df_sets.to_csv(set_txt, sep='\t')#'''

        avg_oob = round(avg_dict['oob'], 1)
        avg_cnt = int(round(avg_dict['count'], 0))
        
        print '\nAverage OOB score: .................... %.1f' % avg_oob
        print '\nAverage number of overlapping sets: ... %s\n' % avg_cnt
        
        print 'Time to make OOB score map: %.1f hours\n' % ((time.time() - t1)/3600)
    
    # Record params in inventory text file
    if 'inventory_txt' in inputs:
        t1 = time.time()
        print 'Getting model info...\n'
        df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
        n_sets = len(df_sets)
        '''if 'sample' in sample_txt:
            n_samples = int(sample_txt.split('_')[1].replace('sample',''))
        inv_columns = df_inv.columns
        if 'n_sets' in inv_columns: df_inv.ix[stamp, 'n_sets'] = n_sets
        if 'n_samples' in inv_columns: df_inv.ix[stamp, 'n_samples'] = n_samples
        if 'support_size' in inv_columns: df_inv.ix[stamp, 'support_size'] = str(support_size)
        if 'sets_per_cell' in inv_columns: df_inv.ix[stamp, 'sets_per_cell'] = sets_per_cell
        if 'max_features' in inv_columns: df_inv.ix[stamp, 'max_features'] = max_features
        info_dir = os.path.dirname(inventory_txt)
        existing_models = fnmatch.filter(os.listdir(info_dir), '%s*' % target_col)
        if len(existing_models) > 0:
            df_inv = df_inv[df_inv.index.isin(existing_models)]#'''


        if 'avg_oob' in inv_columns and make_oob_map: df_inv.ix[stamp, 'avg_oob'] = avg_oob
        if 'avg_count' in inv_columns and make_oob_map: df_inv.ix[stamp, 'avg_count'] = avg_cnt
        if len(df_inv) > 1:
            df_inv.to_csv(inventory_txt, sep='\t')
        else:
            print 'WARNING: Model info not written to inventory_txt...\n' #'''   
        
    
    print 'Total training time: %.1f minutes' % ((time.time() - t0)/60)

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))

#params = '/vol/v2/stem/jdb_test/param_files/train_stem_params.txt'




