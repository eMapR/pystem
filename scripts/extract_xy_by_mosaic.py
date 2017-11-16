'''
Script for extracting data values by (x,y) coordinate from data organized
by TSA.
'''


import os
import fnmatch
import sys
import glob
import shutil
import warnings
#import random
import pandas as pd
import numpy as np
#from gdalconst import *
from scipy import stats
from osgeo import gdal, ogr
import time
from datetime import datetime
from multiprocessing import Pool

import stem
from lthacks import attributes_to_df

gdal.UseExceptions()

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


def find_file(basepath, search_str, tile_str=None, path_filter=None, year=None):
    '''
    Return the full path within the directory tree /baspath/tile_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    if not os.path.exists(basepath.replace('*','')):
        print 'basepath does not exist: \n%s' % basepath
        return None
     
    if tile_str: 
        '''try:
            bp = os.path.join(basepath, tile_str)
        except:
            import pdb; pdb.set_trace()#'''

        # Search the bp directory tree. If search_str is in a file, get the full path.
        paths = []
        for root, dirs, files in os.walk(basepath, followlinks=True):
            if not os.path.basename(root) == tile_str:
                continue
            #these_paths = [os.path.join(root, f) for f in files]
            #these_paths = [f for f in fnmatch.filter(these_paths, '*%s*%s' % (tile_str, search_str)]# if tile_str in f]
            these_paths = [os.path.join(root, f) for f  in fnmatch.filter(files, search_str)]
            paths.extend(these_paths)
    else:
        paths = glob.glob(os.path.join(basepath, search_str))
        
    # If path filter is specified, remove any paths that contain it
    if not path_filter == '':
        [paths.remove(p) for p in paths if fnmatch.fnmatch(p, path_filter)]
        #paths = [p for p in paths if fnmatch.fnmatch(p, path_filter)]
    
    if len(paths) > 1:
        print 'Multiple files found for tile: %s' % tile_str
        #import pdb; pdb.set_trace()
        for p in paths:
            print p
        print 'Selecting the first one found...\n'# '''
        
    if len(paths) < 1:
        import pdb; pdb.set_trace()
        raise IOError(('No files found for with tile {0}, basepath {1}, and ' +\
        'search_str {2}\n').format(tile_str, basepath, search_str))
    
    return paths[0]
    


def within(x, y, poly):
    
    # check if point is a vertex
    if (x,y) in poly: 
        return True

    # Check if point is on a boundary
    n_vertices = len(poly)
    for i in xrange(n_vertices): 
        if i == 0:
            x1, y1 = poly[0]
            x2, y2 = poly[1]
        else:
            x1, y1 = poly[i - 1]
            x2, y2 = poly[i]
        
        # If the sum of the distance from pt1 to pt and pt to pt2 are equal to 
        #   the total distance from pt1 to p2, pt is on the line
        dist1 = ((y1 - y)**2 + (x1 - x)**2)**.5
        dist2 = ((y2 - y)**2 + (x2 - x)**2)**.5
        dist = ((y1 - y2)**2 + (x1 - x2)**2)**.5
        if dist1 + dist2 == dist: # This might not work because of rounding errors
            return True#'''
    '''
    # Only works for square tiles
    xmin = poly[0][0]
    ymax = poly[0][1]
    for i in xrange(1, n_vertices):
        x1, y1 = min(poly[i])
        xmin = min(xmin, x1)
        ymax = max(ymax, x1)
    import pdb; pdb.set_trace()
    if y == ymax or x == xmin:
        return True'''

    '''if (y == y2 and y == y1) or (x == x2 and x == x1):
            return True'''
            
    inside = False
    x1, y1 = poly[0]
    for i in xrange(n_vertices + 1):
        
        x2, y2 = poly[i % n_vertices]
        if y > min(y1, y2): 
            if y <= max(y1, y2):
                if x < max(x1, x2):
                    if y1 != y2:
                        xints = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                        if x1 == x2 or x < xints:
                            inside = ~inside
        x1, y1 = x2, y2

    return inside
        

'''def within(x, y, geometry):
    
    point = ogr.CreateGeometryFromWkt('POINT (%s %s)' % (x, y))
    if geometry.Contains(point): 
        return True
    else:
        return False'''


def par_within(args):
    
    geom_coords, xy_temp, feature_id, n, n_features = args
    pct_progress = float(n + 1)/n_features * 100
    points = [i for i, (x, y) in xy_temp[['x','y']].iterrows() if within(x, y, geom_coords)]
    sys.stdout.write('\rRetreived points for feature %s of %s (%.1f%%)' % (n + 1, n_features, pct_progress))
    sys.stdout.flush()
    return feature_id, points#"""




def extract_from_shp(lyr, df_xy, id_field=None, n_jobs=0):
    ''' Extract the id of the polygon that each point in df_xy is within'''
    
    '''ds = ogr.Open(shp)
    lyr = ds.GetLayer()'''
    
    id_field_given = True
    if not id_field: 
        id_field = 'name'
        id_field_given = False

    if n_jobs:
        t1 = time.time()
        args = []
        n_features = lyr.GetFeatureCount()
        tile_ids = []
        for i, feature in enumerate(lyr):
            #feature = lyr.GetFeature(i)
            geometry = feature.GetGeometryRef()
            geom_coords = stem.get_coords_from_geometry(geometry, multipart='split')
            
            # Initially select only samples that fit within the bounds of the this feature
            min_x, max_x, min_y, max_y = geometry.GetEnvelope()
            xy_temp = df_xy[(df_xy.x >= min_x) &
                            (df_xy.x < max_x) & 
                            (df_xy.y > min_y) & 
                            (df_xy.y <= max_y)] 
            #if id_field_given:
            #    feature_id = feature.GetField(id_field)
            #else:
            #feature_id = feature.GetFID()
            feature_id = feature.GetField(id_field)
            if len(xy_temp) > 0: 
                tile_ids.append(feature_id)
            args.append([geom_coords, xy_temp, feature_id, i, n_features])
            feature.Destroy()
            sys.stdout.write('\rInitial filter of points for (%.1f%%) of features' % (float(i)/n_features * 100))
            sys.stdout.flush()
            #feature.Destroy()
        print '\nTime for filtering: %.1f minutes\n' % ((time.time() - t1)/60)

        # Predict in parallel
        t1 = time.time()
        pool = Pool(n_jobs)
        points = pool.map(par_within, args, 1)
        pool.close()
        pool.join()
        print '\nTime for extraction: %.1f minutes\n' % ((time.time() - t1)/60)
        
        t1 = time.time()
        '''for i, p in points:
            df_xy.ix[p, 'tile_id'] = i
        print 'Time for adding to df: %.1f minutes\n' % ((time.time() - t1)/60)'''
        point_dict = {tile_id: point_ids for tile_id, point_ids in points if len(point_ids) > 0}
    else:
        for feature in lyr:
            #feature = lyr.GetFeature(i)
            geometry = feature.GetGeometryRef()
            
            #geom_coords = stem.get_coords_from_geometry(geometry)
            min_x, max_x, min_y, max_y = geometry.GetEnvelope()
            
            
            # Initially select only samples that fit within the bounds of the this feature
            xy_temp = df_xy[(df_xy.x >= min_x) &
                            (df_xy.x < max_x) & 
                            (df_xy.y > min_y) & 
                            (df_xy.y <= max_y)] 
            
            # With square tiles, this is not necessary
            #points = [i for i, (x, y) in xy_temp[['x','y']].iterrows() if within(x, y, geometry)]
            
            if id_field_given:
                feature_id = feature.GetField(id_field)
            else:
                feature_id = feature.GetFID()
            df_xy.ix[xy_temp.index, 'tile_fid'] = feature_id
            feature.Destroy()
    
    return point_dict
    

def extract_at_xy(df, mosaic, data_tx, val_name):
    
    '''in_cols = df.columns
    x_res, y_res = data_tx[1], data_tx[5]
    ul_x, ul_y = data_tx[0], data_tx[3]
    
    if 'row' not in df.columns:
        df['row'] = [int((y - ul_y)/y_res) for y in df['y']]
        df['col'] = [int((x - ul_x)/x_res) for x in df['x']]'''
    
    n_jobs = 20
    if type(mosaic) == ogr.Layer:
        point_dict  = extract_from_shp(mosaic, df, val_name, n_jobs=n_jobs)
        val_name = 'tile_fid'
     
    # Otherwise, it's a numpy array
    else:
        ar_rows, ar_cols = mosaic.shape
        df[val_name] = [mosaic[r, c] if r > 0 and r < ar_rows and c > 0 and c < ar_cols else None for r,c in zip(df['row'], df['col'])]
        point_dict = {}
    
    return point_dict


def extract_by_rowcol(df, sample_ids, filepath, file_col, var_name, data_band, mosaic_tx, val_cols, data_type, nodata=None, kernel=False):
    '''
    Return a dataframe with values from all pixels defined by row, col in df. 
    If row_dirs and col_dirs are defined, get surrounding pixels in a 
    kernel of size row_dirs ** 1/2, col_dirs ** 1/2.
    '''
    row_dirs = [-1,-1,-1, 0, 0, 0, 1, 1, 1]
    col_dirs = [-1, 0, 1,-1, 0, 1,-1, 0, 1]
    
    #df_file = df[df[file_col] == filepath]
    df_file = df.ix[sample_ids]
    
    # If the filepath is an integer, the file could not be found so just
    #   return bogus values that can be picked out later
    if type(filepath) == int:
        vals = np.full((len(df_file),len(val_cols)), nodata, dtype=np.int32)
        df.ix[sample_ids, var_name] = vals
    
    # Otherwise it should be a real filepath
    else:
        ds = gdal.Open(filepath)
        tx = ds.GetGeoTransform()
        band = ds.GetRasterBand(data_band)
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
       
        #ar = band.ReadAsArray()
        #ds = None
        
        # Calc row, col of the dataset for each point
        row_off = int((tx[3] - mosaic_tx[3])/tx[5])
        col_off = int((tx[0] - mosaic_tx[0])/tx[1])
        df_file = df_file[(df_file.row - row_off < ysize) & (df_file.col - col_off < xsize)]
        
        if kernel:
            if nodata == None:
                nodata = -9999
            
            # Need to buffer the array because of extracting kernel at edges
            ar = np.full((ysize + 2, xsize + 2), nodata, dtype=np.int16)
            ar[1:-1, 1:-1] = band.ReadAsArray()#'''
            row_off -= 1
            col_off -= 1
            
            data_rows = [row - row_off + d for row in df_file.row for d in row_dirs]
            data_cols = [col - col_off + d for col in df_file.col for d in col_dirs]
            
            # Get the data values and stats for each kernel
            try:
                vals = ar[data_rows, data_cols].reshape(len(df_file), len(val_cols))
                
            except:
                import pdb; pdb.set_trace()
        
            # Make a df from the array and from an array of stats for each kernel
            df_vals = pd.DataFrame(vals, columns=val_cols, index=df_file.index)
            df_stats = pd.DataFrame(calc_row_stats(vals, data_type, var_name, nodata), index=df_file.index)

            # Add the new values and stats to the input df 
            df.ix[sample_ids, var_name] = df_stats[var_name]
            stat_cols = list(df_stats.columns)
            df_vals[stat_cols] = df_stats
            df_vals[file_col] = df_file[file_col]
            
        else:
            #import pdb; pdb.set_trace()
            if filepath.endswith('vrt'):
                vals = {}
                for i, (c, r) in df_file[['col', 'row']].iterrows():
                    v = band.ReadAsArray(c - col_off, r - row_off, 1, 1)
                    #import pdb; pdb.set_trace()
                    vals[i] = v[0][0]
                #vals = {var_name: [band.ReadAsArray(c - col_off, r - row_off, 1, 1)[0][0] for i, (c, r) in df_file[['col', 'row']].iterrows()]}
            else:
                ar = band.ReadAsArray()
                try:
                    vals = {var_name: ar[df_file.row - row_off, df_file.col - col_off]}
                except:
                    import pdb; pdb.set_trace()
            df_vals = pd.DataFrame(vals, index=df_file.index)
            try:
                df.ix[df_file.index, var_name] = df_vals[var_name]
            except:
                import pdb; pdb.set_trace()
        band = None
        ds = None
            
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
        val_stats = {var_name: mean, 
                     var_name + '_median': med, 
                     var_name + '_variance': var,
                     var_name + '_stdv': std,
                     var_name + '_range': rnge
                     }
    
    # Otherwise, they're discrete
    else:
        mode, count = stats.mode(ar, axis=1)
        val_stats = {var_name: mode.ravel(), 
                     var_name + '_count': count.ravel()
                     }
        
    return val_stats


def par_extract_by_rowcol(args):
    
    sys.stdout.write('\rFile count: %s' % args[-1])
    sys.stdout.flush()
    return extract_by_rowcol(*args[:-1])
    

def extract_var(year, var_name, by_tile, data_band, data_type, df_tile, df_xy, point_dict, basepath, search_str, path_filter, mosaic_tx, file_count, n_files, nodata=None, kernel=False):
    '''
    Return a dataframe of 
    '''
    t0 = time.time()
    dfs = [] # For storing kernel and stats
    var_col = var_name# + str(year)
    file_col = 'file_' + var_col
    #file_count = last_file
    # Store the filepath for each tile
    if by_tile:
        df_tile[file_col] = [find_file(basepath, search_str.format(year), tile, path_filter)
        for tile in df_tile.tile_id] 
    else:
        df_tile[file_col] = find_file(basepath, search_str.format(year), path_filter=path_filter)
    # Handle any rows for for which the file is null
    if df_tile[file_col].isnull().any():
        df_null = df_tile[df_tile[file_col].isnull()]
        print 'Tiles excluded from extractions for %s from %s:' % (var_name, year)
        for ind, row in df_null.iterrows(): print row['tile_str']
        print ''
        n_null = len(df_null)
        # Make the file name a unique integer so that it can be
        #   distinguished from real files and from other null files
        df_tile.loc[df_null.index, file_col] = range(n_null)
        
    # Get the file string for each xy for each year
    df_xy[file_col] = ''# Creates the column but keeps it empty
    #for tile in df_xy.tile_id.unique():
    for tile in point_dict.keys():
        try:
            #df_xy.loc[df_xy['tile_id'] == tile, file_col] = df_tile.loc[df_tile['tile_id'] == tile, file_col].values[0]
            df_xy.loc[point_dict[tile], file_col] = df_tile.loc[df_tile['tile_id'] == tile, file_col].values[0]
        except:
            import pdb; pdb.set_trace()
    #import pdb; pdb.set_trace()
    # For each file, get the dataset as an array and extract all values at each row col
    val_cols = ['%s_%s' % (var_col, i) for i in range(1,10)]
    #for i, f in enumerate(df_tile[file_col].unique()):
    for i, tile in enumerate(point_dict):
        f = df_tile.ix[df_tile.tile_id == tile, file_col].values[0]
        print 'Extracting for array %s of approximately %s from:\n%s\n'\
        % (file_count, n_files, f)
        dfs.append(extract_by_rowcol(df_xy, point_dict[tile], f, file_col, var_col, data_band,
                                     mosaic_tx, val_cols, data_type, nodata,
                                     kernel
                                     ))
        file_count += 1#'''
        #this_count = i + 1
    '''args = []
    for i, f in enumerate(df_tile[file_col].unique()):
        args.append([df_xy, f, file_col, var_col, data_band, mosaic_tx, val_cols, data_type, nodata, kernel, i + 1 + file_count])
    n_jobs=20
    pool = Pool(n_jobs)
    print 'Extracting from %s-%s files of %s...' % (file_count, file_count + this_count, n_files)
    dfs = pool.map(par_extract_by_rowcol, args, 1)
    pool.close()
    pool.join()'''
    print '\nTime for this variable: %.1f minutes\n' % ((time.time() - t0)/60)
    #file_count += this_count
    
    # Comnbine all the pieces for this year
    df_var = pd.concat(dfs)
    
    return df_var, file_count
    

def main(params, out_dir=None, xy_txt=None, kernel=False, resolution=30, tile_id_field=None):          
    t0 = time.time()
    # Read params. Make variables from each line of the 1-line variables
    inputs, df_vars = read_params(params)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    columns = df_vars.columns
    expected_cols = ['var_name', 'data_type', 'data_band', 'search_str', 'basepath', 'by_tile', 'nodata', 'path_filter']
    #missing_cols = 
    print '\nExtracting variable values for samples:\n%s\n' % xy_txt
    
    if 'out_dir' not in inputs:
        out_dir = os.path.dirname(xy_txt)
    if not os.path.isdir(out_dir):
        #raise IOError('Output directory does not exist: %s' % out_dir)
        print 'Output directory does not exist. Creating new directory at', out_dir
        os.mkdir(out_dir)
    shutil.copy2(params, out_dir)
    
    if not os.path.exists(mosaic_path):
        raise IOError('mosaic_path does not exist: %s' % mosaic_path)
        
    if 'years' in locals(): 
        years =  [int(yr) for yr in years.split(',')]
    else:
        try:
            year_start = int(year_start)
            year_end = int(year_end)
            years = range(year_start, year_end + 1)
        except NameError:
            raise NameError('No list of years or year_start/year_end specified in' +\
            ' param file:\n%s\n. Re-run script with either of these' +\
            ' parameters given.' % params)
    if 'kernel' in inputs:
        if inputs['kernel'].lower().replace('"','') == 'true':
            kernel = True
        else:
            kernel = False
    
    # Get the TSA mosaic as an array
    #if (df_vars.by_tsa == 1).any():
    print 'Reading mosaic dataset...%s\n' % time.ctime(time.time())
    try: 
        mosaic_ds = ogr.Open(mosaic_path)
        if 'resolution' not in inputs:
            warnings.warn('Resolution not specified. Assuming default of 30...\n')
        mosaic = mosaic_ds.GetLayer()
        min_x, max_x, min_y, max_y = mosaic.GetExtent()
        xsize = int((max_x - min_x)/resolution)
        ysize = int((max_y - min_y)/resolution)
        prj = mosaic.GetSpatialRef().ExportToWkt()
        driver = gdal.GetDriverByName('envi')
        x_res = resolution
        y_res = -resolution
        mosaic_tx = min_x, x_res, 0, max_y, 0, y_res
    
    except:
        exc_type, exc_msg, _ = sys.exc_info()
        
        try:
            mosaic_ds = gdal.Open(mosaic_path)
            exc_type_str = str(exc_type).split('.')[-1].replace("'", '').replace('>', '')
            msg = ('Could not open mosaic_path as vector : dataset %s.\n%s: %s' +\
                   '\n\nAttempting to open as raster...\n') % (mosaic_path, exc_type_str, exc_msg)
            warnings.warn(msg)
            
        except:
            exc_type, exc_msg, _ = sys.exc_info()
            exc_type_str = str(exc_type).split('.')[-1].replace("'", '').replace('>', '')
            msg = 'Could not open mosaic_path as vector or raster: dataset %s.\n%s: %s' % (mosaic_path, exc_type_str, exc_msg)
            raise IOError(msg)
        mosaic_tx = mosaic_ds.GetGeoTransform()
        mosaic = mosaic_ds.ReadAsArray()
        #tile_ids = np.unique(df_tile.tile_id) # Assumes there are coords in all TSAs'''


    # Only need to do this once per xy sample set
    # Get the tile id at each xy location and the row and col 
    print 'Extracting tile IDs... %s\n' % time.ctime(time.time())
    df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')
    if np.any(df_vars.by_tile):
        point_dict = extract_at_xy(df_xy, mosaic, mosaic_tx, tile_id_field)# Gets val inplace
        if len(point_dict) > 0:
            tile_ids = point_dict.keys()
        else:
            tile_ids = df_xy[tile_id_field].unique()
    
        #df_mosaic = attributes_to_df(mosaic_path)
        #df_xy = df_xy[~df_xy.tile_fid.isnull()] #drop samples that weren't inside a tile
        #df_xy['tile_id'] = df_mosaic.ix[df_xy['tile_fid'], tile_id_field].tolist()
        #tile_id_field = 'tile_id'
    #import pdb; pdb.set_trace()
    
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = mosaic_tx
    df_xy['row'] = [int((y - ul_y)/y_res) for y in df_xy.y]
    df_xy['col'] = [int((x - ul_x)/x_res) for x in df_xy.x]
    #df_xy = df_xy[df_xy.tile_id > 0]
    #df_xy[tile_id_field] = 0
        
    #tile_ids = df_xy['tile_id'].unique()
    
    if 'tile_txt' in inputs:
        df_tile = pd.read_csv(tile_txt, sep='\t', dtype={'tsa_str':object})
        df_tile['tile_str'] = df_tile['tsa_str']
    else:
        df_tile = pd.DataFrame({'tile_id': tile_ids, 'tile_str': tile_ids})# id string is to ensure compatibility with TSAs

    # For each year, do extractions
    c = 1 #For keeping count of files processed
    n_years = len(years)
    
    #n_vars = len(df_vars)
    #import pdb; pdb.set_trace()
    n_files = (len(df_tile) * len(df_vars[df_vars.by_tile == 1])) * n_years +\
    n_years * (len(df_vars[df_vars.data_band < 0]) + len(df_vars[df_vars.data_band > 0]))
    xy_txt_bn = os.path.basename(xy_txt)
    #out_dir = os.path.join(out_dir, xy_txt_bn.split('.')[0])
    #if not os.path.exists(out_dir): os.mkdir(out_dir)
    
    last_file = 1
    for year in years:
        t1 = time.time()
        df_yr = pd.DataFrame(index=df_xy.index)
        for var_name, var_row in df_vars.iterrows(): #var_name is index col
            # Get variables from row
            search_str  = var_row.search_str
            basepath    = var_row.basepath
            path_filter = var_row.path_filter
            data_type   = var_row.data_type
            by_tile      = var_row.by_tile
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

            df_var, last_file = extract_var(year, var_name, by_tile, data_band,
                                             data_type, df_tile, df_xy, point_dict, basepath,
                                             search_str, path_filter, mosaic_tx,
                                             last_file, n_files, nodata, kernel)
            df_yr[df_var.columns] = df_var
            #last_file += this_count
        
        # Write df_var with all years for this predictor to txt file
        df_xy[df_vars.index.tolist()] = df_yr[df_vars.index.tolist()]
        if kernel:
            this_bn = '%s_%s_kernelstats.txt' % (xy_txt_bn[:-4], year)
            this_txt = os.path.join(out_dir, this_bn)
            df_yr.to_csv(this_txt, sep='\t')

        print 'Time for year %s: %.1f minutes\n\n' % (year, ((time.time() - t1)/60))
        
    mosaic_ds = None
    
    # Write the dataframe to a text file
    #if out_dir:
    #out_cols = [col for col in df_xy.columns if 'file' not in col]
    out_cols = ['x', 'y', target_col] + df_vars.index.tolist()
    
    out_bn = xy_txt_bn.replace('.txt', '_predictors.txt')
    out_txt = os.path.join(out_dir, out_bn)
    #df_out = df_xy.drop(['tile_id', 'row', 'col', 'x', 'y'], axis=1)# May want to keep row/col
    # Remove any rows where anything is null
    null_mask = df_xy.isnull().any(axis=1)
    df_out = df_xy.ix[~null_mask, out_cols]
    df_out.to_csv(out_txt, sep='\t')
    df_xy.ix[null_mask, out_cols].to_csv(out_txt.replace('.txt', '_dropped.txt'), sep='\t')
    print 'Dataframe written to:\n', out_txt       
    print 'Total extraction time: %.1f minutes' % ((time.time() - t0)/60)
    
    return out_txt

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
#df_xy, df_tile = main(params)
    