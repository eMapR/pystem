'''
Script for extracting data values by (x,y) coordinate from data organized
by TSA.
'''

import gdal
import os
import fnmatch
import sys
import glob
#import random
import pandas as pd
import numpy as np
from gdalconst import *
from scipy import stats
import time
from datetime import datetime

# Turn off the stupid setting as copy warning for pandas dataframes
pd.options.mode.chained_assignment = None

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


def find_file(basepath, search_str, tsa_str=None, path_filter=None, year=None):
    '''
    Return the full path within the directory tree /baspath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    if not os.path.exists(basepath.replace('*','')):
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
    
    if len(paths) > 1:
        print 'Multiple files found for tsa: %s' % tsa_str
        import pdb; pdb.set_trace()
        for p in paths:
            print p
        print 'Selecting the first one found...\n'# '''
        
    if len(paths) < 1:
        sys.exit(('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str))
    
    return paths[0]
    
    
"""def find_file(basepath, tsa_str, search_str, path_filter=None):
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
        #these_files = [f for f in files if fnmatch.fnmatch(f, search_str)]
        these_paths = [os.path.join(root, f) for f in files]
        these_paths = fnmatch.filter(these_paths, search_str)
        if len(these_paths) > 0:
            #these_paths = [os.path.join(root, f) for f in these_files]
            paths.extend(these_paths)
    
    # If path filter was specified, remove any paths that contain it
    if not path_filter == '':
        [paths.remove(p) for p in paths if fnmatch.fnmatch(p, path_filter)]
        #paths = [p for p in paths if fnmatch.fnmatch(p, path_filter)]
    
    if len(paths) > 1:
        print 'Multiple files found for tsa: ' + tsa_str
        for p in paths:
            print p
        print 'Selecting the first one found...\n'
        #return None
        
    if len(paths) < 1:
        print ('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str)
        return None
        
    return paths[0]#"""


def extract_at_xy(df, array, data_tx, val_name):
    
    x_res, y_res = data_tx[1], data_tx[5]
    ul_x, ul_y = data_tx[0], data_tx[3]
    ar_rows, ar_cols = array.shape
    
    if 'row' not in df.columns:
        df['row'] = [int((y - ul_y)/y_res) for y in df['y']]
        df['col'] = [int((x - ul_x)/x_res) for x in df['x']]
    
    df[val_name] = [array[r, c] if r > 0 and r < ar_rows and c > 0 and c < ar_cols else None for r,c in zip(df['row'], df['col'])]
    
    return df[['row', 'col', val_name]] 


def extract_by_rowcol(df, filepath, file_col, var_name, data_band, mosaic_tx, val_cols, data_type, nodata=None, kernel=False):
    '''
    Return a df with values from all pixels defined by row, col in df. 
    If row_dirs and col_dirs are defined, get surrounding pixels in a 
    kernel of size row_dirs ** 1/2, col_dirs ** 1/2.
    '''
    row_dirs = [-1,-1,-1, 0, 0, 0, 1, 1, 1]
    col_dirs = [-1, 0, 1,-1, 0, 1,-1, 0, 1]
    
    df_temp = df[df[file_col] == filepath]
    
    # If the filepath is an integer, the file could not be found so just
    #   return bogus values that can be picked out later
    if type(filepath) == int:
        vals = np.full((len(df_temp),len(val_cols)), nodata, dtype=np.int32)
    
    # Otherwise it should be a real filepath
    else:
        ds = gdal.Open(filepath, GA_ReadOnly)
        tx = ds.GetGeoTransform()
        ar = ds.GetRasterBand(data_band).ReadAsArray()
        ds = None
        
        # Calc row, col of the dataset for each point
        row_off = int((tx[3] - mosaic_tx[3])/tx[5])
        col_off = int((tx[0] - mosaic_tx[0])/tx[1])
    
    if kernel:
        data_rows = [row - row_off + d for row in df_temp.row for d in row_dirs]
        data_cols = [col - col_off + d for col in df_temp.col for d in col_dirs]
        
        # Get the data values and stats for each kernel
        vals = ar[data_rows, data_cols].reshape(len(df_temp), len(val_cols))
    
        # Make a df from the array and from an array of stats for each kernel
        df_vals = pd.DataFrame(vals, columns=val_cols, index=df_temp.index)
        df_stats = pd.DataFrame(calc_row_stats(vals, data_type, var_name, nodata), index=df_temp.index)
        
        # Add the new values and stats to the input df 
        df.ix[df[file_col] == filepath, var_name] = df_stats[var_name]
        stat_cols = list(df_stats.columns)
        df_vals[stat_cols] = df_stats
        #df_vals[file_col] = df_temp[file_col]
        
    else:
        vals = {var_name: ar[df_temp.row - row_off, df_temp.col - col_off]}
        df_vals = pd.DataFrame(vals, index=df_temp.index)
        
    return df_vals 


def calc_row_stats(ar, data_type, var_name, nodata):
    '''
    Return a dict of relevant stats for each row in ar given the data_type
    '''
    # Mask out nodata values if the kernel has any
    if nodata != None:
        ar = np.ma.masked_array(ar, ar==nodata, dtype=ar.dtype)
    
    # If the data are continuous
    if data_type == 'continuous':
        mean = np.mean(ar, axis=1, dtype=np.int32)
        med = np.median(ar, axis=1)
        var = np.var(ar, axis=1)
        std = np.std(ar, axis=1)
        rnge = np.ptp(ar, axis=1)
        val_stats = {var_name: mean, var_name + '_median': med, 
                     var_name + '_variance': var, var_name + '_stdv': std,
                     var_name + '_range': rnge
                     }
    
    # Otherwise, they're discrete
    else:
        mode, count = stats.mode(ar, axis=1)
        val_stats = {var_name: mode.ravel(), 
                     var_name + '_count': count.ravel()
                     }
        
    return val_stats

def extract_var(year, var_name, by_tsa, data_band, data_type, df_tsa, df_xy, basepath, search_str, path_filter, mosaic_tx, last_file, n_files, nodata=None):
    '''
    Return a dataframe of 
    '''
    
    dfs = [] # For storing kernel and stats
    var_col = var_name# + str(year)
    file_col = 'file_' + var_col
    file_count = last_file
    # Store the filepath for each TSA
    if by_tsa:
        df_tsa[file_col] = [find_file(basepath, search_str.format(year), tsa, path_filter)
        for tsa in df_tsa.tsa_str] 
    else:
        df_tsa[file_col] = find_file(basepath, search_str.format(year), path_filter=path_filter)
    # Handle any rows for for which the file is null
    if df_tsa[file_col].isnull().any():
        df_null = df_tsa[df_tsa[file_col].isnull()]
        print 'TSAs excluded from extractions for %s from %s:' % (var_name, year)
        for ind, row in df_null.iterrows(): print row['tsa_str']
        print ''
        n_null = len(df_null)
        # Make the file name a unique integer so that it can be
        #   distinguished from real files and from other null files
        df_tsa.loc[df_null.index, file_col] = range(n_null)
        
    # Get the file string for each xy for each year
    df_xy[file_col] = ''# Creates the column but keeps it empty
    for tsa in df_xy.tsa_id.unique():
        df_xy.loc[df_xy['tsa_id'] == tsa, file_col] = df_tsa.loc[df_tsa['tsa_id'] == tsa, file_col].values[0]

    # For each file, get the dataset as an array and extract all values at each row col
    val_cols = ['%s_%s' % (var_col, i) for i in range(1,10)]
    for f in df_tsa[file_col].unique():
        print 'Extracting for array %s of approximately %s from:\n%s\n'\
        % (last_file, n_files, f)
        # If extract from the dataset depending on how year is stored
        dfs.append(extract_by_rowcol(df_xy, f, file_col, var_col, data_band,
                                     mosaic_tx, val_cols, data_type, nodata))
        last_file += 1
        
    # Comnbine all the pieces for this year
    df_var = pd.concat(dfs)
    
    return df_var, last_file - file_count
    

def main(params, out_dir=None, xy_txt=None):          
    t0 = time.time()
    # Read params. Make variables from each line of the 1-line variables
    inputs, df_vars = read_params(params)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    
    print '\nExtracting variable values for samples:\n%s\n' % xy_txt
    
    if not os.path.exists(out_dir):
        print 'Output directory does not exist:\n', out_dir
        return None
    
    if 'years' in locals(): 
        years =  [int(yr) for yr in years.split(',')]
    else:
        try:
            year_start = int(year_start)
            year_end = int(year_end)
            years = range(year_start, year_end + 1)
        except NameError:
            print 'No list of years or year_start/year_end specified in' +\
            ' param file:\n%s\n. Re-run script with either of these' +\
            ' parameters given.' % params
            return None
    
    # Get the TSA mosaic as an array
    #if (df_vars.by_tsa == 1).any():
    print 'Reading mosaic dataset...%s\n' % time.ctime(time.time())
    mosaic_ds = gdal.Open(tsa_mosaic, GA_ReadOnly)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    mosaic_ar = mosaic_ds.ReadAsArray()
        
    
    # Read in the TSA txt and XY txt as dataframes
    print 'Reading in text files... %s\n' % time.ctime(time.time())
    df_tsa = pd.read_csv(tsa_txt, sep='\t', dtype={'tsa_str':object})
    df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')

    # Only need to do this once per xy sample set
    # Get the TSA at each xy location and the row and col 
    print 'Extracting tsa array values... %s\n' % time.ctime(time.time())
    extract_at_xy(df_xy, mosaic_ar, mosaic_tx, 'tsa_id')# Gets val inplace
    df_xy = df_xy[df_xy.tsa_id > 0]
    
    # For each year, do extractions
    tsa_ids = np.unique(df_tsa.tsa_id) # Assumes there are coords in all TSAs
    c = 1 #For keeping count of files processed
    n_years = len(years)
    #n_vars = len(df_vars)
    n_files = (len(df_tsa) * len(df_vars[df_vars.by_tsa == 1])) * n_years +\
    n_years * (len(df_vars[df_vars.data_band < 0]) + len(df_vars[df_vars.data_band > 0]))
    xy_txt_bn = os.path.basename(xy_txt)
    #out_dir = os.path.join(out_dir, xy_txt_bn.split('.')[0])
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    
    last_file = 1
    for year in years:
        t1 = time.time()
        df_yr = pd.DataFrame()
        for var_name, var_row in df_vars.iterrows(): #var_name is index col
            # Get variables from row
            search_str  = var_row.search_str
            basepath    = var_row.basepath
            path_filter = var_row.path_filter
            data_type   = var_row.data_type
            by_tsa      = var_row.by_tsa
            path_filter = ''
            #data_band   = var_row.data_band 
            if np.isnan(var_row.nodata):
                nodata = None
            else:
                nodata = int(var_row.nodata)
            if int(var_row.data_band) < 0: 
                data_band = year - 1984 + 1
            else: 
                data_band = int(var_row.data_band)

            df_var, this_count = extract_var(year, var_name, by_tsa, data_band, data_type, df_tsa, df_xy, basepath, search_str, path_filter, mosaic_tx, last_file, n_files, nodata)
            df_yr[df_var.columns] = df_var
            last_file += this_count
        
        # Write df_var with all years for this predictor to txt file
        this_bn = '%s_%s_kernelstats.txt' % (xy_txt_bn[:-4], year)
        this_txt = os.path.join(out_dir, this_bn)
        df_yr.to_csv(this_txt, sep='\t')
        print 'Time for year %s: %.1f minutes\n\n' % (year, ((time.time() - t1)/60))
        
    mosaic_ds = None
    
    # Write the dataframe to a text file
    #if out_dir:
    out_cols = [col for col in df_xy.columns if 'file' not in col]
    
    out_bn = xy_txt_bn.replace('.txt', '_predictors.txt')
    out_txt = os.path.join(out_dir, out_bn)
    #df_out = df_xy.drop(['tsa_id', 'row', 'col', 'x', 'y'], axis=1)# May want to keep row/col
    # Remove any rows where anything is null
    df_out = df_xy.ix[~df_xy.isnull().any(axis=1), out_cols] 
    df_out.to_csv(out_txt, sep='\t')
    print 'Dataframe written to:\n', out_txt       
    print 'Total extraction time: %.1f minutes' % ((time.time() - t0)/60)

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params)) #'''
    
    
''' testing '''
#tsa_mosaic = '/vol/v1/general_files/datasets/spatial_data/calorewash_TSA_nobuffer.bsq'
#tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'
#xy_txt = '/vol/v2/stem/canopy/samples/canopy_sample3000_20160122_1600.txt'
#basepath = '/vol/v1/proj/lst/outputs/models/randomforest/rfprediction'
#search_strs = {'landcover': 'lst_run1_prediction_voting_lulc_RF_*%s.tif'}
#path_filter = None
#out_dir= '/vol/v2/stem/scripts/testing/'
#data_type = 'discrete'
#years = [2001, 2002]

#df_xy = main(tsa_mosaic, tsa_txt, years, basepath, search_strs, path_filter, 1, data_type, out_dir, xy_txt=xy_txt)
#params = '/vol/v2/stem/param_files/extract_xy_by_tsa_params_topovars.txt'
#df_xy, df_tsa = main(params)
    