'''
Script for extracting data values by (x,y) coordinate from data organized
by TSA.
'''

import gdal
import os
import fnmatch
import random
import pandas as pd
import numpy as np
from gdalconst import *
from scipy import stats
import time
import generate_gsrd.py as gsrd
# Turn off the stupid setting as copy warning for pandas dataframes
pd.options.mode.chained_assignment = None


def find_file(basepath, tsa_str, search_str, path_filter=None):
    '''
    Return the full path within the directory tree /baspath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    if not os.path.exists(basepath):
        print 'basepath does not exist: \n%s' % basepath
        return None
     
    bp = os.path.join(basepath, tsa_str)

    # Search the bp directory tree. If search_str is in a file, get the full path.
    paths = []
    for root, dirs, files in os.walk(bp, followlinks=True):
        these_files = [f for f in files if fnmatch.fnmatch(f, search_str)]
        if len(these_files) > 0:
            these_paths = [os.path.join(root, f) for f in these_files]
            paths.extend(these_paths)
    
    if path_filter:
        paths = [p for p in paths if path_filter not in p]
    
    if len(paths) > 1:
        print 'Multiple files found for tsa: %s\n' % tsa_str
        for p in paths:
            print p
        return None
        
    if len(paths) < 1:
        print ('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str)
        return np.nan
    return paths[0]


def extract_at_xy(df, array, data_tx, nodata):
    
    x_res, y_res = data_tx[1], data_tx[5]
    ul_x, ul_y = data_tx[0], data_tx[3]
    ar_rows, ar_cols = array.shape
    
    if 'row' not in df.columns:
        df['row'] = [int((y - ul_y)/y_res) for y in df['y']]
        df['col'] = [int((x - ul_x)/x_res) for x in df['x']]
    
    df['val'] = [array[r, c] if r > 0 and r < ar_rows and c > 0 and c < ar_cols else None for r,c in zip(df['row'], df['col'])]
    
    return df[['row', 'col', 'val']] 


def extract_by_rowcol(df, filepath, data_band, mosaic_tx, val_cols, data_type, row_dirs=None, col_dirs=None):
    '''
    Return a df with values from all pixels defined by row, col in df. 
    If row_dirs and col_dirs are defined, get surrounding pixels in a 
    kernel of size row_dirs ** 1/2, col_dirs ** 1/2.
    '''
    ds = gdal.Open(filepath, GA_ReadOnly)
    tx = ds.GetGeoTransform()
    ar = ds.GetRasterBand(data_band).ReadAsArray()
    ar_rows, ar_cols = ar.shape
    ds = None
    
    # Calc row, col of the dataset for each point
    row_off = int((tx[3] - mosaic_tx[3])/tx[5])
    col_off = int((tx[0] - mosaic_tx[0])/tx[1])
    data_rows = [row - row_off + d for row in df['row'] for d in row_dirs]
    data_cols = [col - col_off + d for col in df['col'] for d in col_dirs]
    
    # Get the data values and stats for each kernel
    vals = ar[data_rows, data_cols].reshape(len(df), len(val_cols))
    df_vals = pd.DataFrame(vals, columns=val_cols, index=df.index)
    df_stats = pd.DataFrame(calc_row_stats(vals, data_type), index=df.index)
    
    df[val_cols] = df_vals
    df[df_stats.columns] = df_stats
    
    return df


def calc_row_stats(ar, data_type):
    '''
    Return a dict of relevant stats from ar given the data_type
    '''
    # If the data are continuous
    if data_type == 'continuos':
        mean = np.mean(ar, axis=1)
        med = np.median(ar, axis=1)
        var = np.var(ar, axis=1)
        std = np.std(ar, axis=1)
        rnge = np.ptp(ar, axis=1)
        val_stats = {'val': mean, 'median': med, 'variance': var,
                 'stdv': std, 'range': rnge}
    
    # Otherwise, they're discrete
    else:
        mode, count = stats.mode(ar, axis=1)
        val_stats = {'val': mode.ravel(), 'count': count.ravel()}
        
    return val_stats
                

def main(tsa_mosaic, tsa_txt, basepath, search_str, path_filter, data_band, data_type, nodata, out_txt=None, xy_txt=None):      
    
    # Get the TSA mosaic as an array
    print 'Reading mosaic dataset...%s\n' % time.ctime(time.time())
    mosaic_ds = gdal.Open(tsa_mosaic, GA_ReadOnly)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    mosaic_ar = mosaic_ds.ReadAsArray()
    
    # Read in the TSA txt and XY txt as dataframes
    print 'Reading in text files... %s\n' % time.ctime(time.time())
    df_tsa = pd.read_csv(tsa_txt, sep='\t', dtype={'tsa_str':object})
    #df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')
    if xy_txt: 
        df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')
    else:
        #Generate random samples  
    
    # Store the filepath for each TSA
    print 'Finding data files... %s\n' % time.ctime(time.time())
    df_tsa['file'] = [find_file(basepath, tsa, search_str, path_filter) for tsa in df_tsa['tsa_str']]
    #df_tsa['file'] = df_tsa['tsa_str'].apply(lambda x: find_file(basepath, x, search_str, path_filter))
        
    
    # Remove any rows for for which the file is null
    if df_tsa['file'].isnull().any():
        df_null = df_tsa[df_tsa['file'].isnull()]
        print 'TSAs excluded from extractions:'
        for ind, row in df_null.iterrows(): print row['tsa_str']
        print ''
        df_tsa = df_tsa.drop(df_null.index)
    
    # Only need to do this once per xy sample set
    # Get the TSA at each xy location and the row and col 
    print 'Extracting tsa array values... %s\n' % time.ctime(time.time())
    df_xy[['row','col', 'tsa_id']] = extract_at_xy(
    df_xy, mosaic_ar, mosaic_tx, nodata)#'''
    
    # Get the file string for each xy
    print 'Getting file string for each xy... %s\n' % time.ctime(time.time())
    tsa_ids = np.unique(df_tsa['tsa_id'])
    df_xy['file'] = ''
    for tsa in tsa_ids:
        df_xy.loc[df_xy['tsa_id'] == tsa, 'file'] = df_tsa.loc[
        df_tsa['tsa_id'] == tsa, 'file'].values[0]
    print time.ctime(time.time())
    
    # For each file, get the dataset as an array and extract all values at each row col
    print 'Extracting values for each point... %s\n' % time.ctime(time.time())
    row_dirs = [-1,-1,-1, 0, 0, 0, 1, 1, 1]
    col_dirs = [-1, 0, 1,-1, 0, 1,-1, 0, 1]
    val_cols = ['val%s' % i for i in range(1,10)]
    dfs = []
    c = 1
    for f in df_tsa['file']:
        print 'Extracting for file %s of %s:\n%s\n' % (c, len(df_tsa), f)
        df_temp = df_xy[df_xy['file'] == f]
        dfs.append(extract_by_rowcol(df_temp, f, data_band, 
                                        mosaic_tx, val_cols, data_type,
                                        row_dirs, col_dirs))
        c += 1
    df_out = pd.concat(dfs) #Comnbine all the pieces
    mosaic_ds = None
    
    # If out_text is specified, write the dataframe to a text file
    if out_txt:
        df_out.to_csv(out_txt, sep='\t', index=False)
        print 'Dataframe written to:\n', out_txt
    
    return df_out
    
    
    
''' testing '''
tsa_mosaic = '/vol/v1/general_files/datasets/spatial_data/calorewash_TSA_nobuffer.bsq'
tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'
xy_txt = '/vol/v2/stem/scripts/xy_prj.txt'
basepath = '/vol/v1/proj/lst/outputs/models/randomforest/rfprediction'
search_str = 'lst_run1_prediction_voting_lulc_RF_*2002.tif'
path_filter = None
out_txt= '/vol/v2/stem/scripts/testing/extract_xy_test.txt'
data_type = 'discrete'

#df_xy = main(tsa_mosaic, tsa_txt, xy_txt, basepath, search_str, path_filter, 1, data_type, 255, out_txt)
    