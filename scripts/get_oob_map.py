# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 19:54:05 2016

@author: shooper
"""

import time
import warnings
from osgeo import gdal
from gdalconst import *
#from scipy import stats
import pandas as pd
import numpy as np

import stem

gdal.UseExceptions()
warnings.filterwarnings('ignore')


def get_oob_samples(xy_txt, df_train, df_sets):
    '''
    Get the out of bag samples for each set in df_sets. df_train should come from the
    txt file of samples from the model dir. xy_txt should be the txt with 'predictor' 
    at the end of the filename from the sample dir.
    '''
    df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')
    #df_train = pd.read_csv(train_txt, sep='\t')
    #df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    
    # Get tuples of set_id, df of obs within set
    oob_obs = [(i, df_xy[(df_xy['x'] > r[['ul_x', 'lr_x']].min()) &
    (df_xy['x'] < r[['ul_x', 'lr_x']].max()) &
    (df_xy['y'] > r[['ul_y', 'lr_y']].min()) &
    (df_xy['y'] < r[['ul_y', 'lr_y']].max())]) for i, r in df_sets.iterrows()]
    oob_dfs = []
    for i, df in oob_obs:
        df['set_id'] = i
        oob_dfs.append(df) # Get the set id
    df_xy = pd.concat(oob_dfs)
    del oob_obs, oob_dfs        
    
    # Get the out of bag obs by looping trhough each set ID and finding all of
    #   the obs that are within the set and don't share a commmon ID with any
    #   train obs from the set
    oob_dfs = []
    for set_id in df_sets.index:
        this_train = df_train[df_train.set_id == set_id]
        this_xy = df_xy[(df_xy.set_id == set_id)]
        this_xy = this_xy[~this_xy.index.isin(this_train.obs_id)]
        oob_dfs.append(this_xy)
    df_oob = pd.concat(oob_dfs)
    del oob_dfs
    
    return df_oob
    

def main(xy_txt, train_txt, set_txt, err_threshold, val_col, var_cols, ysize, xsize, nodata, mask, n_tiles, tx, support_size, out_dir, file_stamp, prj, driver, inventory_txt=None, oob_txt=None):
    
    df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    
    t0 = time.time()
    if not oob_txt:
        df_train = pd.read_csv(train_txt, sep='\t')
        print '\nGetting OOB samples... ', time.ctime(t0)
        df_oob = get_oob_samples(xy_txt, df_train, df_sets)
        print '%.1f seconds\n' % (time.time() - t0)
    else:
        print '\nReading oob_txt: %s\n' % oob_txt
        df_oob = pd.read_csv(oob_txt, sep='\t')
    
    print 'Getting OOB score map...\n'
    ar_oob, ar_cnt, df_sets = stem.oob_map(ysize, xsize, nodata, mask, n_tiles, tx, support_size, df_oob, df_sets, val_col, var_cols, err_threshold, out_dir, file_stamp, prj, driver)
    df_sets.to_csv(set_txt, sep='\t')
    
    if inventory_txt:
        df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
        avg_oob = np.mean(ar_oob[mask])
        avg_cnt = int(round(np.mean(ar_cnt[mask]), 0))
        df_inv.ix[file_stamp, 'oob_score'] = avg_oob
        df_inv.ix[file_stamp, 'avg_nsets'] = avg_cnt
        df_inv.to_csv(inventory_txt, sep='\t')
        print '\nAverage OOB score: .................... %.1f' % avg_oob
        print '\nAverage number of overlapping sets: ... %s'   % avg_cnt
        
   

   
'''################Testing################'''
"""mosaic_path = '/vol/v1/general_files/datasets/spatial_data/CAORWA_TSA_lt_only.bsq'
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
  
xy_txt = '/vol/v2/stem/imperv//samples/imperv_sample141806_20160702_1714/imperv_sample141806_20160702_1714_predictors.txt'
train_txt = '/vol/v2/stem/imperv/models/imperv_20160928_1216/imperv_20160928_1216_train.txt'
set_txt = '/vol/v2/stem/imperv/models/imperv_20160928_1216/decisiontree_models/imperv_20160928_1216_support_sets.txt'
err_threshold = 10
val_col = 'imperv'
var_cols = ['aspectNESW', 'aspectNWSE', 'brightness', 'delta_bright', 'delta_green', 'delta_nbr', 'delta_wet', 'elevation', 'greenness', 'mse', 'nbr', 'slope', 'time_since', 'wetness']
oob_txt = '/vol/v2/stem/imperv/models/imperv_20160928_1216/imperv_20160928_1216_oob.txt'
inventory_txt = '/vol/v2/stem/imperv/models/model_info.txt'
nodata = 0
n_tiles = 25, 15 
support_size = 400000, 300000
out_dir = os.path.dirname(train_txt)
file_stamp = os.path.basename(out_dir)
 
main(xy_txt, train_txt, set_txt, err_threshold, val_col, var_cols, ysize, xsize, nodata, mask, n_tiles, tx, support_size, out_dir, file_stamp, prj, driver, inventory_txt=inventory_txt, oob_txt=oob_txt)#"""
 
