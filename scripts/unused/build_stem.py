# -*- coding: utf-8 -*-
"""
-generate random xy locations
-sample each predictor with xtract at xy
-generate gsrd and get sample locations within each support set
-for each support set:
    -generate decision tree
    -generate mosaic for each variable
    -predict with decision tree
 

@author: shooper
"""
import os
import sys
import time
import fnmatch
import glob
import random
import shutil
#from itertools import count as itertoolscount
#from random import sample as randomsample
#from string import ascii_lowercase
from osgeo import gdal
from gdalconst import *
from sklearn import tree
#from multiprocessing import Pool
from datetime import datetime
import cPickle as pickle
import pandas as pd
import numpy as np

# Import ancillary scripts
import generate_gsrd as gsrd
import mosaic_by_tsa as mosaic
import aggregate_stem as aggr


"""_data_name_cands = (
    '_data_' + ''.join(randomsample(ascii_lowercase, 10))
    for _ in itertoolscount())

class ForkedData(object):
    '''
    Class used to pass data to child processes in multiprocessing without
    really pickling/unpickling it. Only works on POSIX.

    Intended use:
        - The master process makes the data somehow, and does e.g.
            data = ForkedData(the_value)
        - The master makes sure to keep a reference to the ForkedData object
          until the children are all done with it, since the global reference
          is deleted to avoid memory leaks when the ForkedData object dies.
        - Master process constructs a multiprocessing.Pool *after*
          the ForkedData construction, so that the forked processes
          inherit the new global.
        - Master calls e.g. pool.map with data as an argument.
        - Child gets the real value through data.value, and uses it read-only.
    '''
    
    def __init__(self, val):
        g = globals()
        self.name = next(n for n in _data_name_cands if n not in g)
        g[self.name] = val
        self.master_pid = os.getpid()

    @property
    def value(self):
        return globals()[self.name]

    def __del__(self):
        if os.getpid() == self.master_pid:
            del globals()[self.name]"""


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


def vars_to_numbers(cell_size, support_size, sets_per_cell, min_obs, pct_train, n_tiles, nodata):
    '''
    Return variables as ints or floats
    '''
    cell_size = [int(i) for i in cell_size.split(',')]
    support_size = [int(i) for i in support_size.split(',')]
    sets_per_cell = int(sets_per_cell)
    min_obs = int(min_obs)
    pct_train = float(pct_train)
    n_tiles = [int(i) for i in n_tiles.split(',')]
    nodata = int(nodata)
    
    return cell_size, support_size, sets_per_cell, min_obs, pct_train, n_tiles, nodata


def find_file(basepath, search_str, tsa_str=None, path_filter=None):
    '''
    Return the full path within the directory tree /baspath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    if not os.path.exists(basepath):
        print 'basepath does not exist: \n%s' % basepath
        return None
     
    if tsa_str: 
        bp = os.path.join(basepath, tsa_str)

        # Search the bp directory tree. If search_str is in a file, get the full path.
        paths = []
        for root, dirs, files in os.walk(bp, followlinks=True):
            these_paths = [os.path.join(root, f) for f in files]
            these_paths = fnmatch.filter(these_paths, search_str)
            paths.extend(these_paths)
    else:
        paths = glob.glob(os.path.join(basepath, search_str))
        
    # If path filter is specified, remove any paths that contain it
    if not path_filter == '':
        [paths.remove(p) for p in paths if fnmatch.fnmatch(p, path_filter)]
        #paths = [p for p in paths if fnmatch.fnmatch(p, path_filter)]
    
    '''if len(paths) > 1:
        print 'Multiple files found for tsa: ' + tsa_str
        for p in paths:
            print p
        print 'Selecting the first one found...\n'# '''
        
    if len(paths) < 1:
        sys.exit(('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str))
    
    return paths[0]


def fit_tree(x_train, y_train):
    ''' '''
    dt = tree.DecisionTreeClassifier()
    dt.fit(x_train, y_train)
    
    return dt


def write_decisiontree(dt, filename):
    ''' 
    Pickle a decision tree and write it to filename. Return the filename.
    '''
    with open(filename, 'w+') as f:
        pickle.dump(dt, f)
        
    return filename
    

def write_model(out_dir, df_sets):
    '''
    Write STEM decision trees and dataframe of locations to disk
    '''
    this_dir = os.path.join(out_dir, 'decisiontree_models')
    if not os.path.exists(this_dir):
        os.mkdir(this_dir)
    
    stamp = os.path.split(out_dir)[1]
    dt_bn = stamp + '_decisiontree_%s'
    
    dt_file = os.path.join(this_dir, dt_bn)
    df_sets['dt_file'] = [write_decisiontree(row.dt_model, dt_file % set_id)\
                         for set_id, row in df_sets.iterrows()]
    
    set_txt = os.path.join(this_dir, stamp + '_support_sets.txt')
    df_sets['set_id'] = df_sets.index
    df_sets.drop('dt_model', axis=1).to_csv(set_txt, sep='\t', index=False)
    
    print 'Support set dataframe and decision trees written to:\n', this_dir
    

def get_predict_array(args):

    ar = mosaic.get_mosaic(*args[1:])
    return args[0], ar.ravel()
    
    
def get_predictors(df_var, mosaic_tx, tsa_strs, tsa_ar, ar_coords, save_stuff=None):
    '''
    Return an array of flattened predictor arrays where each predictor is a 
    separate column
    '''
    # #### tsa_ar shoud be a keyword set to None
    t0 = time.time()
    predictors = []
    #rows, cols = tsa_ar.shape
    #df_predict = np.empty((rows * cols, len(df_var)), dtype=np.int32)
    for ind, var in enumerate(df_var.index):
        this_tsa_ar = np.copy(tsa_ar)
        data_band, search_str, basepath, by_tsa, path_filter = df_var.ix[var]
        if by_tsa:
            files = [find_file(basepath, search_str, tsa, path_filter) for tsa in tsa_strs]
            ar_var = mosaic.get_mosaic(mosaic_tx, tsa_strs, this_tsa_ar, ar_coords, data_band, files)
            ''' delete this stuff'''
            '''m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
            ul_x, ul_y = ar_coords[['ul_x', 'ul_y']]
            tx = ul_x, 30.0, 0.0, ul_y, 0, -30.0
            set_id, mosaic_dir, prj, driver = save_stuff
            out_path = mosaic_dir + '/%s_%s.bsq' % (var, set_id)
            mosaic.array_to_raster(ar_var, tx, prj, driver, out_path, GDT_Int32)#'''
            
        else:
            this_file = find_file(basepath, search_str, path_filter=path_filter)
            tx_, ar_var, roff_, coff_ = mosaic.get_array(this_file, data_band, ar_coords)
        #predictors[var] = ar_var.ravel()
        predictors.append(ar_var.ravel())

    #df_predict = pd.DataFrame(predictors)
    ar = np.vstack(predictors).T

    del predictors
    
    print 'Finished getting arrays...', time.time() - t0        
    return ar
    
    
def predict_set(set_id, df_var, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, dt, nodata, save_stuff=None):
    '''
    Return a predicted array for set, set_ind
    '''
    # Get an array of tsa_ids within the bounds of ar_coords
    tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, 1, ar_coords, mosaic_tx,
                                            xsize, ysize, nodata=nodata)
    tsa_ar[tsa_ar==0] = nodata
    # Get the ids of TSAs this kernel covers
    tsa_ids = np.unique(tsa_ar)
    tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=nodata]
    array_shape = tsa_ar.shape

    # Get an array of predictors where each column is a flattened 2D array of a
    #   single predictor variable
    ar_predict = get_predictors(df_var, mosaic_tx, tsa_strs, tsa_ar, ar_coords)#, save_stuff)
    del tsa_ar #Release resources from the tsa array
    
    t0 = time.time()
    nodata_mask = np.any(~(ar_predict==nodata), axis=1)
    predictions = dt.predict(ar_predict[nodata_mask]).astype(np.int32)
    ar_prediction = np.full(ar_predict.shape[0], nodata, dtype=np.int32)
    ar_prediction[nodata_mask] = predictions
    
    print 'Finished predicting...', time.time() - t0
    
    return ar_prediction.reshape(array_shape)


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
    df_overlap = df_sets[
    ((df_sets.ctr_x - tile_bounds.ctr_x).abs() < max_x_dist) &
    ((df_sets.ctr_y - tile_bounds.ctr_y).abs() < max_y_dist)]
    
    return df_overlap


def calc_offset_from_tile(tile_ul, array_ul, tx):
    '''
    Return the row and col offset of a data array from a tsa_array
    '''
    tile_x, tile_y = tile_ul
    array_x, array_y = array_ul
    
    row_off = int((array_y - tile_y)/tx[5])
    col_off = int((array_x - tile_x)/tx[1])
    
    #return pd.Series((row_off, col_off))
    return row_off, col_off


def fill_tile_band(tile_size, tile_coords, set_coords, ar_pred, tx, nodata):
    '''
    Fill an array of zeros of shape tile_size, located at tile_coords with an 
    offset array, ar_pred, located at set_coords
    '''
    # Calc offsets
    row_off, col_off = calc_offset_from_tile(tile_coords[['ul_x','ul_y']],
                                             set_coords[['ul_x','ul_y']],
                                             tx)
    # Get the offset indices of each array 
    tile_inds, set_inds = mosaic.get_offset_array_indices(
        tile_size, ar_pred.shape, (row_off, col_off))
    tile_row_u, tile_row_d, tile_col_l, tile_col_r = tile_inds
    set_row_u,  set_row_d,  set_col_l,  set_col_r  = set_inds
    
    # Fill just the part of the array that overlaps
    ar_tile = np.full(tile_size, np.nan)
    ar_pred = ar_pred.astype(float)
    ar_pred[ar_pred == nodata] = np.nan
    try:
        ar_tile[tile_row_u:tile_row_d, tile_col_l:tile_col_r] =\
        ar_pred[set_row_u:set_row_d, set_col_l:set_col_r]
    except Exception as e:
        print e
        print '\nProblem with offsets'
        print row_off, col_off, set_coords, tile_coords      

    return ar_tile


def main(params):
    
    '''### copy params to out_dir #### '''
    
    #read_params(params)
    inputs, df_var = read_params(params)

    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    try:
        num_vars = vars_to_numbers(cell_size, support_size, sets_per_cell,
                                   min_obs, pct_train, n_tiles, nodata)
        cell_size, support_size, sets_per_cell, min_obs, pct_train, n_tiles, nodata = num_vars
        str_check = sample_txt, target_col, mosaic_path, tsa_txt, dep_var_name, out_dir
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
        return None
    
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    stamp = '{0}_{1}_{2}'.format(dep_var_name, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it exists
    shutil.copy2(params, out_dir) #Copy the params for reference
    
    # Get samples and support set bounds
    if 'gsrd_shp' not in locals(): gsrd_shp = None
    out_txt = os.path.join(out_dir, stamp + '.txt')
    dfs = gsrd.get_gsrd(mosaic_path, cell_size, support_size, sets_per_cell,
                        sample_txt, min_obs, pct_train, dep_var_name, out_txt,
                        gsrd_shp)
    df_train, df_test, df_sets = dfs
    support_sets = df_train.set_id.unique()

    # Check that df_train has exactly the same columns as variables specified in df_vars
    #   Last four characters in each column of df_train should be year
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
    df_sets['dt_model'] = [fit_tree(x_train.ix[x_train.set_id==s, predict_cols],\
    y_train.ix[y_train.set_id==s, target_col]) for s in support_sets]
    
    # Write df_sets and each decison tree to disk
    write_model(out_dir, df_sets)
    
    mosaic_ds = gdal.Open(mosaic_path, GA_ReadOnly)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    xsize = mosaic_ds.RasterXSize
    ysize = mosaic_ds.RasterYSize
    prj = mosaic_ds.GetProjection()
    driver = mosaic_ds.GetDriver()
    
    t0 = time.time()
    
    predict_dir = os.path.join(out_dir, 'predctions')
    os.mkdir(predict_dir)
    # Loop through each set and generate predictions
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    c = 1
    total_sets = len(support_sets)
    predictions = {}
    for set_id, row in df_sets.iterrows():
        print 'Predicting for set %s of %s' % (c, total_sets)
        ar_coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
        ar_predict = predict_set(set_id, df_var, mosaic_ds, ar_coords, 
                                 mosaic_tx, xsize, ysize, row.dt_model, nodata)
        #predictions[set_id] = ar_predict
        
        tx = ar_coords['ul_x'], x_res, x_rot, ar_coords['ul_y'], y_rot, y_res
        out_path = predict_dir + '/prediction_%s.bsq' % set_id
        mosaic.array_to_raster(ar_predict, tx, prj, driver, out_path, GDT_Int32, nodata=nodata)
        c += 1
    mosaic_ds = None                  
    print '\nTotal time for predictions: %.1f minutes' % ((time.time() - t0)/60)#'''
    
    #Aggregate predictions by tile and stitch them back together
    aggr.aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_tx, support_size, predict_dir, df_sets, out_dir, stamp, prj, driver)


def predict_set_from_disk(df_sets, set_id, params):
    
    inputs, df_var = read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    df_var = df_var.reindex(df_var.index.sort_values())
    this_set = df_sets.ix[set_id]
    with open(this_set.dt_file, 'rb') as f: 
        dt_model = pickle.load(f)
    
    mosaic_ds = gdal.Open(mosaic_path, GA_ReadOnly)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    xsize = mosaic_ds.RasterXSize
    ysize = mosaic_ds.RasterYSize
    prj = mosaic_ds.GetProjection()
    driver = mosaic_ds.GetDriver()
    
    ar_coords = this_set[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
    mosaic_dir = '/vol/v2/stem/canopy/canopy_20160212_2016/var_mosaics'
    saving_stuff = set_id, mosaic_dir, prj, driver
    ar_predict = predict_set(set_id, df_var, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, dt_model, saving_stuff)
    return ar_predict
    
    '''out_dir = '/vol/v2/stem/scripts/testing'
    out_path = os.path.join(out_dir, 'predict_rerun_%s.bsq' % set_id)
    
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    tx = this_set.ul_x, x_res, x_rot, this_set.ul_y, y_rot, y_res
    mosaic.array_to_raster(ar_predict, tx, prj, driver, out_path, GDT_Int32)'''

''' ############# Testing ################ '''
#sample_txt = '/vol/v2/stem/canopy/samples/canopy_sample3000_20160122_1600_predictors.txt'
#target_col = 'value'
#mosaic_path = '/vol/v1/general_files/datasets/spatial_data/CAORWA_TSA_lt_only.bsq'
#tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'
#cell_size = (300000, 200000)
#support_size = (400000, 300000)
#sets_per_cell = 10
#min_obs = 25
#pct_train = .63
#dep_var_name = 'canopy'
#n_tiles = 10
#out_dir = '/vol/v2/stem/canopy/models/'
#set_id, ar, df_sets, df_train = main(sample_txt, target_col, mosaic_path, cell_size, support_size, sets_per_cell, min_obs, pct_train, dep_var_name, n_tiles, out_dir)

params = '/vol/v2/stem/param_files/build_stem_params_nomse.txt'
#predictions, df_sets, df_train = main(params)
stuff = main(params)
'''set_txt = '/vol/v2/stem/canopy/outputs/canopy_20160212_2016/decisiontree_models/canopy_20160212_2016_support_sets.txt'
df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
tsa_ar = predict_set_from_disk(df_sets, 341, params)'''

'''tile_size = [size * 30 for size in tile_size]
shp = '/vol/v2/stem/extent_shp/orwaca.shp'
coords, extent = gsrd.get_coords(shp)
gsrd.plot_sets_on_shp(coords, 900, sets[20][1], (400000, 300000), df_tiles.ix[sets[20][0]], tile_size)'''
#for s in sets:
#    print 'Number of sets: ', len(s)
#    coords, extent = gsrd.get_coords(shp)
#    gsrd.plot_sets_on_shp(coords, 500, s, support_size)
