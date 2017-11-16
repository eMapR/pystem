# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:22:37 2017

@author: shooper
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import cPickle as pickle
from osgeo import gdal
from glob import glob

import stem
import mosaic_by_tsa as mosaic
from lthacks import array_to_raster

def main(model_dir, n_tiles, **kwargs):
    
    t0 = time.time()
    
    n_tiles = [int(n) for n in n_tiles.split(',')]
    if not os.path.isdir(model_dir):
        message = 'model directory given does not exist or is not a directory: ', model_dir
        raise IOError(message)
    
    model = os.path.basename(model_dir)
    dt_dir = os.path.join(model_dir, 'decisiontree_models')
    set_txt = os.path.join(dt_dir, '%s_support_sets.txt' % model)
    df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')

    pred_param_path = glob(os.path.join(model_dir, 'predict_stem_*params.txt'))[0]
    predict_params, df_var = stem.read_params(pred_param_path)
    train_param_path = glob(os.path.join(model_dir, 'train_stem_*params.txt'))[0]
    train_params, _ = stem.read_params(train_param_path)
    df_var.sort_index(inplace=True)
    
    nodata = int(predict_params['nodata'].replace('"', ''))   
    if len(kwargs) == 0:
        var_ids = df_sets.max_importance.unique()
        var_names = df_var.ix[var_ids].index
        variables = zip(var_ids, var_names)
    else:
        variables = [(variable_id, variable_name) for variable_name, variable_id in kwargs]
    
        
    mask_path = os.path.join(model_dir, '%s_vote.bsq' % model)
    if not os.path.exists(mask_path):
        mask_path = mask_path.replace('.bsq', '.tif')
    mask_ds = gdal.Open(mask_path)
    mask_tx = mask_ds.GetGeoTransform()
    xsize = mask_ds.RasterXSize
    ysize = mask_ds.RasterYSize
    prj = mask_ds.GetProjection()
    df_tiles, df_tiles_rc, tile_size = stem.get_tiles(n_tiles, xsize, ysize, mask_tx)
    total_tiles = len(df_tiles)
    df_tiles['tile'] = df_tiles.index
    
    # Find the tiles that have only nodata values
    t1 = time.time()
    print '\nFinding empty tiles...'
    mask = mask_ds.ReadAsArray() == nodata
    empty_tiles = stem.find_empty_tiles(df_tiles, ~mask, mask_tx)
    mask_ds = None
    print '%s empty tiles found of %s total tiles\n%.1f minutes\n' %\
    (len(empty_tiles), total_tiles, (time.time() - t1)/60)
    # Select only tiles that are not empty
    df_tiles = df_tiles.select(lambda x: x not in empty_tiles)
    total_tiles = len(df_tiles)
    
    #some_set = df_sets.iloc[0]
    support_size = [int(s) for s in train_params['support_size'].replace('"','').split(',')]
    set_size = [int(abs(s/mask_tx[1])) for s in support_size]
    
    out_dir = os.path.join(model_dir, 'importance_maps')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    print variables
    for vi, (v_id, v_name) in enumerate(variables):

        t1 = time.time()
        print 'Making map for %s: %s of %s variables\n' % (v_name, vi + 1, len(variables))

        ar = np.full((ysize, xsize), nodata, dtype=np.uint8)
        
        for i, (t_ind, t_row) in enumerate(df_tiles.iterrows()):
            t2 = time.time()
            print 'Aggregating for %s of %s tiles' % (i + 1, total_tiles)
            
            # Calculate the size of this tile in case it's at the edge where the
            #   tile size will be slightly different
            this_size = abs(t_row.lr_y - t_row.ul_y), abs(t_row.lr_x - t_row.ul_x)
            df_these_sets = stem.get_overlapping_sets(df_sets, t_row, this_size, support_size)
            
            rc = df_tiles_rc.ix[t_ind]
            this_size = rc.lr_r - rc.ul_r, rc.lr_c - rc.ul_c
            n_sets = len(df_these_sets)
            
            # Load overlapping predictions from disk and read them as arrays
            tile_ul = t_row[['ul_x','ul_y']]
            
            print n_sets, ' Overlapping sets'
            importance_bands = []
            
            importance_values = []
            for s_ind, s_row in df_these_sets.iterrows():
                
                # Calculate offset and array/tile indices
                offset = stem.calc_offset(tile_ul, (s_row.ul_x, s_row.ul_y), mask_tx)
                #if abs(offset[0]) > this_size[0] or abs(offset[1] > this_size[1]):
                
                tile_inds, a_inds = mosaic.get_offset_array_indices(tile_size, set_size, offset)
                
                # Get feature with maximum importance and fill tile with that val
                try:
                    with open(s_row.dt_file, 'rb') as f: 
                        dt_model = pickle.load(f)
                    importance_value = int(dt_model.feature_importances_[v_id] * 100)
                    importance_values.append(importance_value)
                    #filled = np.full((nrows, ncols), importance_value, dtype=np.uint8)
                    #import_band = stem.fill_tile_band(this_size, filled, tile_inds, nodata)
                    import_band = np.full(this_size, np.nan, dtype=np.float16)
                    import_band[tile_inds[0]:tile_inds[1], tile_inds[2]:tile_inds[3]] = importance_value
                    importance_bands.append(import_band)
                except Exception as e:
                    print e
                    continue#'''
            
            print 'Average importance for this tile: %.1f' % np.mean(importance_values)
            #Aggregate
            importance_stack = np.dstack(importance_bands)
            importance_tile = np.nanmean(importance_stack, axis=2)
            tile_mask = mask[rc.ul_r : rc.lr_r, rc.ul_c : rc.lr_c] | np.isnan(importance_tile)
            importance_tile[tile_mask] = nodata
            ar[rc.ul_r : rc.lr_r, rc.ul_c : rc.lr_c] = np.round(importance_tile).astype(np.uint8)
            print 'Aggregation time for this tile: %.1f minutes\n' % ((time.time() - t2)/60)
            
            '''temp_dir = os.path.join(out_dir, 'delete')
            if not os.path.isdir(temp_dir):
                os.mkdir(temp_dir)
            t_tx = tile_ul[0], 30, 0, tile_ul[1], 0, -30
            array_to_raster(np.round(importance_tile).astype(np.uint8), t_tx, prj, gdal.GetDriverByName('gtiff'), os.path.join(temp_dir, 'delete_%s.tif' % t_ind), gdal.GDT_Byte, 255, True)'''
        out_path = os.path.join(out_dir, '%s_importance_%s.tif' % (model, v_name))
        try:
            array_to_raster(ar, mask_tx, prj, gdal.GetDriverByName('gtiff'), out_path, gdal.GDT_Byte, nodata)
        except Exception as e:
            print e
            import pdb; pdb.set_trace()
        print 'Time for this variable: %.1f minutes\n' % ((time.time() - t1)/60)
    
    print '\nTotal time for %s variables: %.1f hours\n' % (len(variables), ((time.time() - t0)/3600))


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
        
        
        
            
            
            
    

    
    
    
    
    
    
    