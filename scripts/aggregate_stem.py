# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 20:46:47 2016

@author: shooper
"""

import os
import sys
import time
import warnings
import traceback
from multiprocessing import Pool
from osgeo import gdal
from gdalconst import *
#from scipy import stats
import pandas as pd
import numpy as np
import cPickle as pickle

import mosaic_by_tsa as mosaic

gdal.UseExceptions()
warnings.filterwarnings('ignore')

def read_params(txt):
    '''
    Return a dictionary from parsed parameters in txt
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
    n_skip_lines = 0 #Keep track of the number of lines w/o a ";"
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))
            n_skip_lines +=1
    
    # Get the lines with information about each variable as a df
    skip_lines = range(len(input_vars) - n_skip_lines, len(input_vars))
    df_vars = pd.read_csv(txt, sep='\t', index_col='var_name', skip_blank_lines=True, skiprows=skip_lines)
    # Drop any rows for which basepath or search str are empty
    df_vars.dropna(inplace=True, subset=['basepath','search_str'])
    df_vars.fillna({'path_filter': ''}, inplace=True)
    
    print 'Parameters read from:\n', txt, '\n'
    return d, df_vars
    
    
def get_tiles(n_tiles, xsize, ysize, tx=None):
    '''
    Return a dataframe representing a grid defined by bounding coords.
    Tiles have rows and cols defined by n_tiles and projected coords 
    defined by tx.
    '''
    tile_rows = ysize/n_tiles[0]
    tile_cols = xsize/n_tiles[1]
    
    # Calc coords by rows and columns
    ul_rows = np.tile([i * tile_rows for i in range(n_tiles[0])], n_tiles[1]) 
    ul_cols = np.repeat([i * tile_cols for i in range(n_tiles[1])], n_tiles[0])
    lr_rows = ul_rows + tile_rows
    lr_cols = ul_cols + tile_cols
    # Make sure the last row/col lines up with the dataset 
    lr_rows[-1] = ysize
    lr_cols[-1] = xsize
    ctr_rows = ul_rows + tile_rows/2
    ctr_cols = ul_cols + tile_cols/2
    
    coords = {'ul_c': ul_cols, 'ul_r': ul_rows,
              'lr_c': lr_cols, 'lr_r': lr_rows,
              }
    df_rc = pd.DataFrame(coords).reindex(columns=['ul_r', 'lr_r', 'ul_c', 'lr_c'])
    
    #if tx: #If the coords need to be projected, not returned as row/col
    # Calc projected coords 
    ul_x = ul_cols * tx[1] + tx[0]
    ul_y = ul_rows * tx[5] + tx[3]
    lr_x = lr_cols * tx[1] + tx[0]
    lr_y = lr_rows * tx[5] + tx[3]
    ctr_x = ctr_cols * tx[1] + tx[0]
    ctr_y = ctr_rows * tx[5] + tx[3]
    
    coords_prj = {'ul_x': ul_x, 'ul_y': ul_y,
                  'lr_x': lr_x, 'lr_y': lr_y,
                  'ctr_x': ctr_x, 'ctr_y': ctr_y
                  }
    
    df_prj = pd.DataFrame(coords_prj, dtype=int)
    df_prj = df_prj.reindex(columns=['ul_x', 'ul_y', 'lr_x', 'lr_y', 'ctr_x', 'ctr_y'])
    
    return df_prj, df_rc, (tile_rows, tile_cols)


def get_overlapping_sets(df_sets, tile_bounds, tile_size, support_size):
    '''
    Return a dataframe of support sets that overlap the tile defined by
    tile bounds
    '''
    set_ysize, set_xsize = support_size
    tile_ysize, tile_xsize = tile_size
    
    # Calculate the max distance that the centers of the tile and the set
    #   could be from one another
    max_x_dist = set_xsize/2 + tile_xsize/2
    max_y_dist = set_ysize/2 + tile_ysize/2
    
    # Get sets where both the x and y center coords are within the respective
    #   max distances of one another
    overlap = df_sets[
    ((df_sets.ctr_x - tile_bounds.ctr_x).abs() < max_x_dist) &
    ((df_sets.ctr_y - tile_bounds.ctr_y).abs() < max_y_dist)]
    #import pdb; pdb.set_trace()
    
    return overlap


def calc_offset(tile_ul, array_ul, tx):
    '''
    Return the row and col offset of a data array from a tsa_array
    '''
    tile_x, tile_y = tile_ul
    array_x, array_y = array_ul
    
    row_off = int((array_y - tile_y)/tx[5])
    col_off = int((array_x - tile_x)/tx[1])
    
    #return pd.Series((row_off, col_off))
    return row_off, col_off


def load_predictions(p_dir, df_sets, tile_ul, tile_size):

    predictions = {}
    
    for set_id in df_sets.index:
        f = os.path.join(p_dir, 'prediction_%s.bsq' % set_id)
        ds = gdal.Open(f)
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        tx = ds.GetGeoTransform()
        offset = calc_offset(tile_ul, (tx[0], tx[3]), tx)
        t_inds, a_inds = mosaic.get_offset_array_indices(tile_size, (ysize, xsize), offset)
        nrows = a_inds[1] - a_inds[0]
        ncols = a_inds[3] - a_inds[2]
        ar = ds.ReadAsArray(a_inds[2], a_inds[0], ncols, nrows)
        #if nrows < 1 or ncols < 1:
        #    import pdb; pdb.set_trace()
        predictions[set_id] = ar, t_inds
        ds = None
    
    return predictions
    

def fill_tile_band(tile_size, ar_pred, tile_inds, nodata):
    '''
    Fill an array of zeros of shape tile_size, located at tile_coords with an 
    offset array, ar_pred, located at set_coords
    '''
    # Fill just the part of the array that overlaps
    try:
        ar_tile = np.full(tile_size, np.nan)
        ar_pred = ar_pred.astype(float)
        ar_pred[ar_pred == nodata] = np.nan

        ar_tile[tile_inds[0]:tile_inds[1], tile_inds[2]:tile_inds[3]] = ar_pred
        #ar_pred[set_row_u:set_row_d, set_col_l:set_col_r]
    except Exception as e:
        import pdb; pdb.set_trace()
        print e
        print '\nProblem with offsets'
        print tile_inds      

    return ar_tile


def get_max_importance(dt):
    
    importance = dt.feature_importances_
    ind = np.argmax(importance)
  
    return ind, importance


def important_features(dt, ar, nodata):
    ''' 
    Return an array of size ar.size where each pixel is the feature from dt 
    that is most commonly the feature of maximum importance for the assigned 
    class of that pixel
    '''
    
    features = dt.tree_.feature
    mask = features >= 0 #For non-nodes, feature is arbitrary
    features = features[mask]
    feature_vals = np.unique(features)
    values = dt.tree_.value
    # Mask out non-nodes and reshape. 1th dimension is always 1 for single
    #   classification problems
    values = values[mask, :, :].reshape(len(features), values.shape[2])
    
    # Loop through each feature and get count of leaf nodes for each class
    sum_list = []
    for f in feature_vals:
        these_sums = np.sum(values[features == f, :], axis=0)
        sum_list.append(these_sums)
    sums = np.vstack(sum_list)
    feat_inds = np.argmax(sums, axis=0)
    classes = dt.classes_
    max_features = {c: f for c,f in zip(classes, feat_inds)}

    # Map the features to the data values
    max_val = np.max(classes)
    mp = np.arange(0, max_val + 1)
    mp[classes] = [max_features[c] for c in classes]
    t3 = time.time()
    ar_out = np.full(ar.shape, nodata, dtype=np.int32)
    t4 = time.time()
    data_mask = ar != nodata
    ar_out[data_mask] = mp[ar[data_mask]]
    
    return ar_out
  
  
def find_empty_tiles(df, nodata_mask, tx):
    
    empty = []
    for i, row in df.iterrows():
        ul_r, ul_c = calc_offset((tx[0], tx[3]), row[['ul_x', 'ul_y']], tx)
        lr_r, lr_c = calc_offset((tx[0], tx[3]), row[['lr_x', 'lr_y']], tx)
        this_mask = nodata_mask[ul_r : lr_r, ul_c : lr_c]
        
        if not this_mask.any():
            empty.append(i)
    
    return empty


def mode(ar, axis=0, nodata=-9999):
    ''' 
    Code from internet to get mode along given axis faster than stats.mode()
    '''
    if ar.size == 1:
        return (ar[0],1)
    elif ar.size == 0:
        raise Exception('Attempted to find mode on an empty array!')
    try:
        axis = [i for i in range(ar.ndim)][axis]
    except IndexError:
        raise Exception('Axis %i out of range for array with %i dimension(s)' % (axis,ar.ndim))

    srt = np.sort(ar, axis=axis)
    dif = np.diff(srt, axis=axis)
    shape = [i for i in dif.shape]
    shape[axis] += 2
    indices = np.indices(shape)[axis]
    index = tuple([slice(None) if i != axis else slice(1,-1) for i in range(dif.ndim)])
    indices[index][dif == 0] = 0
    indices.sort(axis=axis)
    bins = np.diff(indices, axis=axis)
    location = np.argmax(bins, axis=axis)
    mesh = np.indices(bins.shape)
    index = tuple([slice(None) if i != axis else 0 for i in range(dif.ndim)])
    index = [mesh[i][index].ravel() if i != axis else location.ravel() for i in range(bins.ndim)]
    counts = bins[tuple(index)].reshape(location.shape)
    index[axis] = indices[tuple(index)]
    modals = srt[tuple(index)].reshape(location.shape)
    
    return modals#, counts


def weighted_mean(ar, b, c=5, a=1):
    '''
    Calculate the Gaussian weighted mean of a 3D array. Gaussian curve equation: 
    f(x) = ae ** -((x - b)**2/(2c ** 2)), where a adjusts the height of the 
    curve, b adjusts position along the x axis, and c adjusts the width (stdv)
    of the curve.
    '''
    try:
        b = np.dstack([b for i in range(ar.shape[-1])])
        gaussian = (a * np.e) ** -((np.float_(ar) - b)**2/(2 * c ** 2))
        sums_2d = np.nansum(gaussian, axis=(len(ar.shape) - 1))
        sums = np.dstack([sums_2d for i in range(ar.shape[-1])])
        weights = gaussian/sums
        w_mean = np.nansum(ar * weights, axis=(len(ar.shape) - 1))
    except:
        sys.exit(traceback.print_exception(*sys.exc_info()))
    
    return np.round(w_mean,0).astype(np.int32)
    

def get_nonforest_mask(lc_path, ag_path, lc_vals, ag_vals):
    '''
    Create a boolean mask where lc_path==lc_vals or ag_path==ag_vals.
    '''
    # Read in the datasets as arrays
    ds_lc = gdal.Open(lc_path)
    tx_lc = ds_lc.GetGeoTransform()
    ar_lc = ds_lc.GetRasterBand(1).ReadAsArray()
    ds_lc = None
    
    ds_ag = gdal.Open(ag_path)
    tx_ag = ds_ag.GetGeoTransform()
    ar_ag = ds_ag.ReadAsArray()
    ds_ag = None   
    
    # Calc offsets
    offset = calc_offset((tx_lc[0], tx_lc[1]), (tx_ag[0], tx_ag[1]), tx_ag)
    lc_inds, ag_inds = mosaic.get_offset_array_indices(ar_lc.shape, ar_ag.shape, offset)
    
    # In case these were read from a text file, integerize them
    try:
        lc_vals = np.array([int(v) for v in lc_vals.split(',')])
        ag_vals = np.array([int(v) for v in ag_vals.split(',')])
    except:
        pass
    
    # Get masks and combine them
    #lc_mask = np.in1d(ar_lc[lc_inds[0]:lc_inds[1], lc_inds[2]:lc_inds[3]], lc_vals)
    mask = np.in1d(ar_lc.ravel(), lc_vals).reshape(ar_lc.shape)
    lc_view = mask[lc_inds[0]:lc_inds[1], lc_inds[2]:lc_inds[3]]
    ag_mask = np.in1d(ar_ag[ag_inds[0]:ag_inds[1], ag_inds[2]:ag_inds[3]], ag_vals)
    ag_shape = ag_inds[1] - ag_inds[0], ag_inds[3] - ag_inds[2]
    mask[lc_inds[0]:lc_inds[1], lc_inds[2]:lc_inds[3]] = np.logical_or(lc_view, ag_mask.reshape(ag_shape))
    
    del ar_lc, ar_ag, ag_mask
    
    return mask, tx_lc


def mask_array(ar, mask, tx_ar, tx_mask, mask_val=0):
    '''
    Set ar == 0 where mask is true
    '''
    
    ar_ul = tx_ar[0], tx_ar[3]
    mask_ul = tx_mask[0], tx_mask[3]
    offset = calc_offset(ar_ul, mask_ul, tx_ar)
    
    a_inds, m_inds = mosaic.get_offset_array_indices(ar.shape, mask.shape, offset)
    view = ar[a_inds[0]:a_inds[1], a_inds[2]:a_inds[3]]
    view[mask[m_inds[0]:m_inds[1], m_inds[2]:m_inds[3]]] = mask_val


def par_mean(ar):
    ''' Helper function for parallelizing mean '''
    return np.nanmean(ar, axis=(len(ar.shape) - 1))


def par_mode(ar):
    ''' Helper function for parallelizing mode '''
    return mode(ar, axis=(len(ar.shape) - 1))


def par_stdv(ar):
    ''' Helper function for parallelizing standard deviation '''
    return np.nanstd(ar, axis=(len(ar.shape) - 1))


def par_sum(ar):
    ''' Helper function for parallelizing sum '''
    return np.sum(ar, axis=(len(ar.shape) - 1))

def par_wmean(args):
    ''' Helper function for parallelizing weighted_mean '''
    ar, vote, c = args
    #try:
    return weighted_mean(ar, vote, c=c)
    #except:
        #traceback.print_exception(*sys.exc_info())
    

"""def aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_tx, support_size, prediction_dir, df_sets, out_dir, file_stamp, prj, driver):
    
    t0 = time.time()
    ar_mean = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_vote = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_stdv = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_coun = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_impr = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_wtmn_10 = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_wtmn_20 = np.full((ysize, xsize), nodata, dtype=np.int32)#'''
    
    df_tiles, df_tiles_rc, tile_size = get_tiles(n_tiles, xsize, ysize, mosaic_tx)
    total_tiles = n_tiles[0] * n_tiles[1]
    
    
    for t_ind, t_row in df_tiles.ix[:1, :].iterrows():
        t1 = time.time()
        print 'Aggregating for %s of %s tiles' % (t_ind + 1, total_tiles)
        
        # Calculate the size of this tile in case it's at the edge where the
        #   tile size will be slightly different
        rc = df_tiles_rc.ix[t_ind]
        this_size = rc.lr_r - rc.ul_r, rc.lr_c - rc.ul_c
        
        df_these_sets = get_overlapping_sets(df_sets, t_row, this_size, support_size)        
        n_sets = len(df_these_sets)
        set_ids = df_these_sets.index.tolist()
        
        # Load overlapping predictions from disk and read them as arrays
        #files = [os.path.join(prediction_dir, 'prediction_%s.bsq' % set_id)  for set_id in set_ids]
        #predictions = load_predictions(files)
        tile_ul = t_row[['ul_x','ul_y']]
        predictions = load_predictions(prediction_dir, df_these_sets, tile_ul, this_size)
        
        # If there are no overalpping sets for this tile, fill the tile with nodata
        if n_sets == 0:
            print 'No overlapping sets for this tile'
            this_mean = np.full(this_size, nodata, dtype=np.int32)
            this_vote = np.full(this_size, nodata, dtype=np.int32)
            this_stdv = np.full(this_size, nodata, dtype=np.int32)
            this_coun = np.full(this_size, nodata, dtype=np.int32)
            this_impr = np.full(this_size, nodata, dtype=np.int32)
            this_wtmn_10 = np.full(this_size, nodata, dtype=np.int32)
            this_wtmn_20 = np.full(this_size, nodata, dtype=np.int32)
        
        # Otherwise, aggregate all overlapping sets for each pixel
        else:
            print n_sets, ' Overlapping sets'
            t2 = time.time()
            pred_bands = []
            importance_bands = []
            for s_ind, s_row in df_these_sets.iterrows():
                s_row = df_these_sets.ix[s_ind]
                
                # Fill tile with prediction
                ar_pred, tile_inds = predictions[s_ind]
                pred_band = fill_tile_band(this_size, ar_pred, tile_inds, nodata)
                pred_bands.append(pred_band)
                
                # Get feature with maximum importance and fill tile with that val
                with open(s_row.dt_file, 'rb') as f: 
                    dt_model = pickle.load(f)
                max_importance = get_max_importance(dt_model)
                ar_import = np.full(ar_pred.shape, max_import, dtype=np.int32)
                import_band = fill_tile_band(this_size, ar_import, tile_inds, nodata)
                importance_bands.append(import_band)
            print 'Filling tiles: %.1f seconds' % ((time.time() - t2))
            
            ar_tile = np.dstack(pred_bands)
            nd_impr = np.dstack(importance_bands)
            del pred_bands, importance_bands, ar_import
            t3 = time.time()
            n_workers = 40
            p = Pool(n_workers)
            chunksize = ar_tile.shape[0]/n_workers
            this_mean = np.vstack(p.map(par_mean, ar_tile, chunksize))
            this_vote = np.vstack(p.map(par_mode, ar_tile, chunksize))
            this_stdv = np.vstack(p.map(par_stdv, ar_tile, chunksize))
            this_coun = np.vstack(p.map(par_sum, ~np.isnan(ar_tile), chunksize))
            this_impr = np.vstack(p.map(par_mode, nd_impr, chunksize))
            this_wtmn_10 = np.vstack(p.map(par_wmean, [(a, v, 10) for a, v in zip(np.array_split(ar_tile, n_workers), np.array_split(this_vote, n_workers))], chunksize))
            this_wtmn_20 = np.vstack(p.map(par_wmean, [(a, v, 20) for a, v in zip(np.array_split(ar_tile, n_workers), np.array_split(this_vote, n_workers))], chunksize))
            p.close()
            
            print this_mean.shape
            nans = np.isnan(this_mean)
            this_mean[nans] = nodata
            this_stdv[nans] = nodata
            this_impr[nans] = nodata
            this_vote[nans] = nodata
            this_coun[nans] = nodata
            this_wtmn_10[nans] = nodata
            this_wtmn_20[nans] = nodata
            
            print 'Aggregating: %.1f minutes' % ((time.time() - t3)/60)
        print 'Total aggregation time: %.1f minutes\n' % ((time.time() - t1)/60)    
        
        # Fill in the tile in the final arrays
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[t_ind]
        ar_mean[ul_r : lr_r, ul_c : lr_c] = this_mean.astype(np.int32)
        ar_vote[ul_r : lr_r, ul_c : lr_c] = this_vote.astype(np.int32)
        ar_stdv[ul_r : lr_r, ul_c : lr_c] = this_stdv.astype(np.int32)
        ar_coun[ul_r : lr_r, ul_c : lr_c] = this_coun
        ar_impr[ul_r : lr_r, ul_c : lr_c] = this_impr
        ar_wtmn_10[ul_r : lr_r, ul_c : lr_c] = this_wtmn_10
        ar_wtmn_20[ul_r : lr_r, ul_c : lr_c] = this_wtmn_20
    
    # Mask arrays
    #mask, tx_mask = get_nonforest_mask(lc_path, ag_path, lc_vals, ag_vals)
    #mask_array(ar_mean, mask, mosaic_tx, tx_mask)
    #mask_array(ar_vote, mask, mosaic_tx, tx_mask)
    
    # Write final rasters to disk
    out_path = os.path.join(out_dir, file_stamp + '_mean.bsq')
    mosaic.array_to_raster(ar_mean, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = out_path.replace('mean', 'vote')
    mosaic.array_to_raster(ar_vote, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)   
    
    out_path = out_path.replace('vote', 'stdv')
    mosaic.array_to_raster(ar_stdv, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = out_path.replace('stdv', 'count')
    mosaic.array_to_raster(ar_coun, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = out_path.replace('count', 'importance')
    mosaic.array_to_raster(ar_impr, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)  
    
    out_path = os.path.join(out_dir, file_stamp + '_weightedmean_10.bsq')
    mosaic.array_to_raster(ar_wtmn_10, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = os.path.join(out_dir, file_stamp + '_weightedmean_20.bsq')
    mosaic.array_to_raster(ar_wtmn_20, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)#'''
    
    print '\nTotal aggregation run time: %.1f minutes' % ((time.time() - t0)/60)
    #return predictions, df_sets, df_train
    #return df_tiles, df_these_sets, ar_mean, ar_tile, ar_out
    #return ar_out
    del ar_mean, ar_coun, ar_impr, ar_pred, ar_stdv, ar_tile, ar_vote, ar_wtmn_10, ar_wtmn_20, this_coun, this_impr, this_mean, this_stdv, this_vote, this_wtmn_10, this_wtmn_20"""
    

def aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_path, support_size, prediction_dir, df_sets, out_dir, file_stamp, prj, driver, mosaic_nodata=0):
    
    t0 = time.time()
    ar_mean = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_vote = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_stdv = np.full((ysize, xsize), nodata, dtype=np.int32)
    #ar_coun = np.full((ysize, xsize), nodata, dtype=np.int32)
    ar_impr = np.full((ysize, xsize), nodata, dtype=np.int32)
    #ar_wtmn_10 = np.full((ysize, xsize), nodata, dtype=np.int32)#'''
    #ar_wtmn_20 = np.full((ysize, xsize), nodata, dtype=np.int32)
    
    mosaic_ds = gdal.Open(mosaic_path)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    df_tiles, df_tiles_rc, tile_size = get_tiles(n_tiles, xsize, ysize, mosaic_tx)
    total_tiles = len(df_tiles)
    df_tiles['tile'] = df_tiles.index
    
    
    # Find the tiles that have only nodata values
    t1 = time.time()
    print '\nFinding empty tiles...'
    mask = mosaic_ds.ReadAsArray() != mosaic_nodata
    empty_tiles = find_empty_tiles(df_tiles, mask, mosaic_tx)
    mosaic_ds = None
    mask = None
    print '%s empty tiles found of %s total tiles\n%.1f minutes\n' %\
    (len(empty_tiles), total_tiles, (time.time() - t1)/60)
    # Select only tiles that are not empty
    df_tiles = df_tiles.select(lambda x: x not in empty_tiles)
    total_tiles = len(df_tiles)

    # Get feature importances and max importance per set
    t1 = time.time()
    print 'Getting importance values...'
    importance_list = []
    df_sets['max_importance'] = nodata
    for s, row in df_sets.iterrows():
        with open(row.dt_file, 'rb') as f: 
            dt_model = pickle.load(f)
        max_importance, this_importance = get_max_importance(dt_model)
        df_sets.ix[s, 'max_importance'] = max_importance
        importance_list.append(this_importance)
    importance = np.array(importance_list).mean(axis=0)
    pct_import = importance / importance.sum()
    print '%.1f minutes\n' % ((time.time() - t1)/60)#'''
    
    # For each tile, find overlapping sets and calc mode and/or mean for all 
    #   overlapping sets
    #del_dir = '/vol/v2/stem/imperv/models/delete/overlapping'
    #out_txt = os.path.join(del_dir, 'overlapping_%s.txt')
    for i, (t_ind, t_row) in enumerate(df_tiles.iterrows()):
        t1 = time.time()
        print 'Aggregating for %s of %s tiles' % (i + 1, total_tiles)
        
        # Calculate the size of this tile in case it's at the edge where the
        #   tile size will be slightly different
        this_size = abs(t_row.lr_y - t_row.ul_y), abs(t_row.lr_x - t_row.ul_x)
        df_these_sets = get_overlapping_sets(df_sets, t_row, this_size, support_size)
         
        ''' delete '''
        #df_these_sets.to_csv(out_txt % t_ind, sep='\t')
        #continue
        
        rc = df_tiles_rc.ix[t_ind]
        this_size = rc.lr_r - rc.ul_r, rc.lr_c - rc.ul_c
        n_sets = len(df_these_sets)
        
        # Load overlapping predictions from disk and read them as arrays
        tile_ul = t_row[['ul_x','ul_y']]
        predictions = load_predictions(prediction_dir, df_these_sets, tile_ul, this_size)
        
        # If there are no overalpping sets for this tile, fill the tile with nodata
        if t_ind in empty_tiles:
            print 'No overlapping sets for this tile'
            continue
            this_mean = np.full(this_size, nodata, dtype=np.int32)
            this_vote = np.full(this_size, nodata, dtype=np.int32)
            this_stdv = np.full(this_size, nodata, dtype=np.int32)
            #this_coun = np.full(this_size, nodata, dtype=np.int32)
            this_impr = np.full(this_size, nodata, dtype=np.int32)
            #this_wtmn_10 = np.full(this_size, nodata, dtype=np.int32)
            #this_wtmn_20 = np.full(this_size, nodata, dtype=np.int32)
        
        # Otherwise, aggregate all overlapping sets for each pixel
        else:
            print n_sets, ' Overlapping sets'
            t2 = time.time()
            pred_bands = []
            importance_bands = []
            for s_ind, s_row in df_these_sets.iterrows():
                s_row = df_these_sets.ix[s_ind]
                
                # Fill tile with prediction
                ar_pred, tile_inds = predictions[s_ind]
                pred_band = fill_tile_band(this_size, ar_pred, tile_inds, nodata)
                pred_bands.append(pred_band)
                
                # Get feature with maximum importance and fill tile with that val
                try:
                    with open(s_row.dt_file, 'rb') as f: 
                        dt_model = pickle.load(f)
                    #max_importance, importance = get_max_importance(dt_model)
                    #importance_list.append(importance)
                    ar_import = np.full(ar_pred.shape, s_row.max_importance, dtype=np.int32)
                    #max_importance = important_features(dt_model, ar_pred, nodata)
                    import_band = fill_tile_band(this_size, ar_import, tile_inds, nodata)
                    importance_bands.append(import_band)
                except Exception as e:
                    print e
                    continue
            print 'Filling tiles: %.1f seconds' % ((time.time() - t2))
            
            ar_tile = np.dstack(pred_bands)
            nd_impr = np.dstack(importance_bands)
            del pred_bands, importance_bands#, ar_import
            t3 = time.time()
            this_mean = np.nanmean(ar_tile, axis=2)
            this_vote = mode(ar_tile, axis=2)
            this_stdv = np.nanstd(ar_tile, axis=2) * 100 #Multiply b/c converting to int
            #this_coun = np.sum(~np.isnan(ar_tile), axis=2)
            this_impr = mode(nd_impr, axis=2)
            #this_wtmn_10 = weighted_mean(ar_tile, this_vote, c=10)
            #this_wtmn_20 = weighted_mean(ar_tile, this_vote, c=20)
            
            nans = np.isnan(this_vote)
            this_mean[nans] = nodata
            this_stdv[nans] = nodata
            this_impr[nans] = nodata
            this_vote[nans] = nodata
            #this_coun[nans] = nodata
            #this_wtmn_10[nans] = nodata
            #this_wtmn_20[nans] = nodata
            
            print 'Aggregating: %.1f minutes' % ((time.time() - t3)/60)
        print 'Total aggregation time: %.1f minutes\n' % ((time.time() - t1)/60)    
        
        # Fill in the tile in the final arrays
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[t_ind]
        ar_mean[ul_r : lr_r, ul_c : lr_c] = this_mean.astype(np.int32)
        ar_vote[ul_r : lr_r, ul_c : lr_c] = this_vote.astype(np.int32)
        ar_stdv[ul_r : lr_r, ul_c : lr_c] = this_stdv.astype(np.int32)
        #ar_coun[ul_r : lr_r, ul_c : lr_c] = this_coun
        ar_impr[ul_r : lr_r, ul_c : lr_c] = this_impr
        #ar_wtmn_10[ul_r : lr_r, ul_c : lr_c] = this_wtmn_10
        #ar_wtmn_20[ul_r : lr_r, ul_c : lr_c] = this_wtmn_20
    
    # Mask arrays
    #mask, tx_mask = get_nonforest_mask(lc_path, ag_path, lc_vals, ag_vals)
    #mask_array(ar_mean, mask, mosaic_tx, tx_mask)
    #mask_array(ar_vote, mask, mosaic_tx, tx_mask)
    
    # Write final rasters to disk
    out_template = os.path.join(out_dir, file_stamp + '_%s.bsq')
    out_path = out_template % 'mean'
    mosaic.array_to_raster(ar_mean, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = out_template % 'vote'
    mosaic.array_to_raster(ar_vote, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)   
    
    out_path = out_template % 'stdv'
    mosaic.array_to_raster(ar_stdv, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    #out_path = out_path.replace('stdv', 'countagg')
    #mosaic.array_to_raster(ar_coun, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = out_template % 'importance'
    mosaic.array_to_raster(ar_impr, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)  
    
    #out_path = os.path.join(out_dir, file_stamp + '_weightedmean_10.bsq')
    #mosaic.array_to_raster(ar_wtmn_10, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = os.path.join(out_dir, file_stamp + '_weightedmean_20.bsq')
    #mosaic.array_to_raster(ar_wtmn_20, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    print '\nTotal aggregation run time: %.1f hours' % ((time.time() - t0)/3600)
    #return predictions, df_sets, df_train
    #return df_tiles, df_these_sets, ar_mean, ar_tile, ar_out
    #return ar_out
    del ar_impr, ar_pred, ar_stdv, ar_tile, this_impr, this_mean, this_stdv, this_vote, #this_wtmn_20, ar_coun, this_coun, this_wtmn_10, ar_wtmn_20, ar_wtmn_10"""
    
    return ar_mean, ar_vote, pct_import, df_sets


"""mosaic_path = '/vol/v1/general_files/datasets/spatial_data/CAORWA_TSA_lt_only.bsq'
ds = gdal.Open(mosaic_path)
xsize = ds.RasterXSize
ysize = ds.RasterYSize
mosaic_tx = ds.GetGeoTransform()
prj = ds.GetProjection()
driver = ds.GetDriver()
ds = None

params = '/vol/v2/stem/imperv/models/imperv_20160725_2002/predict_stem_params.txt'
inputs, df_var = read_params(params)
for i in inputs:
    exec ("{0} = str({1})").format(i, inputs[i])
df_var = df_var.reindex(df_var.index.sort_values())

nodata = -9999
support_size = 400000, 300000
n_tiles = 25, 15
prediction_dir = '/vol/v2/stem/imperv/models/imperv_20160725_2002/decisiontree_predictions'
set_txt = '/vol/v2/stem/imperv/models/imperv_20160725_2002/decisiontree_models/imperv_20160725_2002_support_sets.txt'
df_sets = pd.read_csv(set_txt, sep='\t', index_col = 'set_id')
out_dir = os.path.dirname(prediction_dir)
file_stamp = 'imperv_20160725_2002'

aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_tx, support_size, prediction_dir, df_sets, out_dir, file_stamp, prj, driver)#"""

'''lc_path = '/vol/v1/proj/lst/outputs/models/randomforest/rfprediction_mosaic/yearly/lst_run1_prediction_voting_lulc_RF_mosaic_2001.bsq'
ag_path = '/vol/v2/stem/canopy/truth_map/crop_mask/crop_CAORWA_2010/crop_mask_no_trees_2010.tif'
lc_vals = [11, 12]
ag_vals = [1]

pred_path = '/vol/v2/stem/canopy/outputs/canopy_20160217_2202/canopy_20160217_2202_final_vote.bsq'
ds_p = gdal.Open(pred_path)
ar_p = ds_p.ReadAsArray()
tx = ds_p.GetGeoTransform()
prj = ds_p.GetProjection()
driver = ds_p.GetDriver()
out_path = pred_path.replace('.bsq', '_masked.bsq')

mask, tx_mask = get_nonforest_mask(lc_path, ag_path, lc_vals, ag_vals)
mask_array(ar_p, mask, tx, tx_mask)
'''

'''pred_path = '/vol/v2/stem/canopy/outputs/canopy_20160424_0916/predctions/prediction_387.bsq'
set_txt = '/vol/v2/stem/canopy/outputs/canopy_20160424_0916/decisiontree_models/canopy_20160424_0916_support_sets.txt'
#df_sets = pd.read_csv(set_txt, sep='\t')
dt_file = '/vol/v2/stem/canopy/outputs/canopy_20160424_0916/decisiontree_models/canopy_20160424_0916_decisiontree_387'
with open(dt_file, 'rb') as f:
    dt = pickle.load(f)

ds_p = gdal.Open(pred_path)
ar_p = ds_p.ReadAsArray()
tx = ds_p.GetGeoTransform()
prj = ds_p.GetProjection()
driver = ds_p.GetDriver()
out_path = '/vol/v2/stem/scripts/testing/importance.bsq'
nodata = -9999

t0 = time.time()
ar = important_features(dt, ar_p, nodata)
print time.time() - t0
mosaic.array_to_raster(ar, tx, prj, driver, out_path, gdal.GDT_Int32, nodata)#'''
