#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 

@author: shooper
"""
import os
import sys
import time
import fnmatch
import glob
import shutil
import random
import warnings
import traceback
from itertools import count as itertoolscount
from random import sample as randomsample
from string import ascii_lowercase
from osgeo import gdal, ogr, gdal_array
from sklearn import tree, metrics
from multiprocessing import Pool
import matplotlib.pyplot as plt
import cPickle as pickle
import pandas as pd
import numpy as np
import pdb

# Import ancillary scripts
import mosaic_by_tsa
import evaluation as ev
from lthacks import get_min_numpy_dtype, array_to_raster, df_to_shp #using this makes no longer portable to any other system than Islay

gdal.UseExceptions()
warnings.filterwarnings('ignore')
# Turn off the annoying setting value on a copy warning
pd.options.mode.chained_assignment = None

_data_name_cands = (
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
            del globals()[self.name]


def read_params(txt):
    '''
    Return a dictionary and a dataframe from parsed parameters in txt
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
    n_skip_lines = 0 #Keep track of the number of lines w/ a ";"
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))
            n_skip_lines +=1
    
    # Get the lines with information about each variable as a df
    ''' Consider making skip_lines a list of indices rather than a range (maybe search \t)'''
    skip_lines = range(len(input_vars) - n_skip_lines, len(input_vars))
    df_vars = pd.read_csv(txt, sep='\t', index_col='var_name', skip_blank_lines=True, skiprows=skip_lines)
    # Drop any rows for which basepath or search str are empty
    df_vars.dropna(inplace=True, subset=['basepath','search_str'])
    df_vars.fillna({'path_filter': ''}, inplace=True)
    
    print '\nParameters read from:\n', txt, '\n'
    return d, df_vars


def vars_to_numbers(cell_size, support_size, sets_per_cell, min_obs, max_features, pct_train):
    '''
    Return variables as ints or floats
    '''
    cell_size = [int(i) for i in cell_size.split(',')]
    support_size = [int(i) for i in support_size.split(',')]
    sets_per_cell = int(sets_per_cell)
    min_obs = int(min_obs)
    if max_features:
        if '.' in max_features: max_features = float(max_features)
        else:
            try: max_features = int(max_features)
            except: pass
    if pct_train: pct_train = float(pct_train)    
    return cell_size, support_size, sets_per_cell, min_obs, max_features, pct_train


def get_raster_bounds(raster):
    '''
    Return the xy bounds of raster
    '''
    ds = gdal.Open(raster)
    tx = ds.GetGeoTransform()
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    x_size = ds.RasterXSize
    y_size = ds.RasterYSize
    ds = None
    
    lr_x = ul_x + x_size * x_res
    lr_y = ul_y + y_size * y_res
    
    min_x = min(ul_x, lr_x)
    min_y = min(ul_y, lr_y)
    max_x = max(ul_x, lr_x)
    max_y = max(ul_y, lr_y)
    
    return min_x, min_y, max_x, max_y, x_res, y_res, tx
    

def generate_gsrd_grid(cell_size, min_x, min_y, max_x, max_y, x_res, y_res):
    ''' 
    Return a dataframe of bounding coordinates 
    '''
    y_size, x_size = cell_size
    
    # Get a randomly defined coordinate within the study area for a seed
    #  upper left coord
    seed_x = random.randint(min_x, max_x)
    seed_y = random.randint(min_y, max_y)
    

    # Calculate how many grid cells there are on either side of the seed
    #   for both x and y dimensions. +1 or not only works for aea projection
    ncells_less_x = int((seed_x - min_x)/x_size + 1)
    ncells_more_x = int((max_x - seed_x)/x_size)
    ncells_less_y = int((seed_y - min_y)/y_size)
    ncells_more_y = int((max_y - seed_y)/y_size + 1)
    
    # Calculate the ul coordinate of each cell
    ul_x = sorted([seed_x - (i * x_size) for i in range(ncells_less_x + 1)])
    ul_x.extend([seed_x + (i * x_size) for i in range(1, ncells_more_x + 1)])
    ul_y = sorted([seed_y - (i * y_size) for i in range(ncells_less_y + 1)])
    ul_y.extend([seed_y + (i * y_size) for i in range(1, ncells_more_y + 1)])
    
    # Make a list of lists where each sub-list is bounding coords of a cell
    x_res_sign = int(x_res/abs(x_res))
    y_res_sign = int(y_res/abs(y_res))
    cells = [[x,
              y,
              x + (x_size * x_res_sign),\
              y + (y_size * y_res_sign)] for x in ul_x for y in ul_y]

    return cells


def sample_gsrd_cell(n, cell_bounds, x_size, y_size, x_res, y_res, tx, snap_coord=None):
    '''
    Return a list of bounding coordinates for n support sets from 
    randomly sampled x,y center coords within bounds
    '''
    ul_x, ul_y, lr_x, lr_y = cell_bounds
    min_x, max_x = min(ul_x, lr_x), max(ul_x, lr_x)
    min_y, max_y = min(ul_y, lr_y), max(ul_y, lr_y)
    
    # Calculate support set centers and snap them to the ultimate raster grid
    #   Use the snap coordinate if given
    if snap_coord:
        x_remain = snap_coord[0] % x_res
        y_remain = snap_coord[1] % y_res
    # Otherwise just use the ul corner of the dataset bounding box
    else:
        x_remain = (tx[0] % x_res)
        y_remain = (tx[3] % y_res)
    x_centers = [int(x/x_res) * x_res - x_remain for x in random.sample(xrange(min_x, max_x + 1), n)]
    y_centers = [int(y/y_res) * y_res - y_remain for y in random.sample(xrange(min_y, max_y + 1), n)]
    
    # Calculate bounding coords from support set centers and make sure 
    #   they're still snapped
    x_res_sign = int(x_res/abs(x_res))
    y_res_sign = int(y_res/abs(y_res))  
    x_remain = (x_size/2) % x_res
    y_remain = (y_size/2) % y_res
    ul_x_ls = [int(round(x - ((x_size/2 + x_remain) * x_res_sign), 0)) for x in x_centers]
    lr_x_ls = [int(round(x + ((x_size/2 - x_remain) * x_res_sign), 0)) for x in x_centers]
    ul_y_ls = [int(round(y - ((y_size/2 + y_remain) * y_res_sign), 0)) for y in y_centers]
    lr_y_ls = [int(round(y + ((y_size/2 - y_remain) * y_res_sign), 0)) for y in y_centers]

    # Store coords for each support set in a dataframe
    these_bounds = zip(ul_x_ls, ul_y_ls, lr_x_ls, lr_y_ls, x_centers, y_centers)
    columns = ['ul_x', 'ul_y', 'lr_x', 'lr_y', 'ctr_x', 'ctr_y']
    df = pd.DataFrame(these_bounds, columns=columns)
    
    return df


def split_train_test(df, pct_train=.8):
    '''
    Return two dataframes: one of training samples and the other testing. 
    '''
    # Get unique combinations of row/col
    unique = [(rc[0],rc[1]) for rc in df[['row','col']].drop_duplicates().values]
    n_train = int(len(unique) * pct_train)
    
    # Randomly sample n_training locations
    try:
        rowcol = random.sample(unique, n_train) #Sample unique row,col
        rows = [rc[0] for rc in rowcol]
        cols = [rc[1] for rc in rowcol]
        # Get all of the obs from those row cols for training points
        df_train = df[df.row.isin(rows) & df.col.isin(cols)]
        # Get all of the obs not from those row,cols for testing
        df_test = df[~df.index.isin(df_train.index)]

    except ValueError:
        # A value error likely means n_train > len(df). This isn't
        #   really possible with the current code, but whatever.
        return None, None
    
    return df_train, df_test
 

def get_obs_within_sets(df_train, df_sets, min_obs, pct_train=None):
    '''
    Return dfs of training and testing obs containing all points within
    each set in df_sets. For training, only return sets if the they 
    contain >= min_obs.
    '''
    # Split train and test sets if pct_train is specified
    df_test = pd.DataFrame()
    if pct_train:
        df_train, df_test = split_train_test(df_train, pct_train)
    
    # Get a list of (index, df) tuples where df contains only points
    #   within the support set    
    obs = [(i, df_train[(df_train['x'] > r[['ul_x', 'lr_x']].min()) &
    (df_train['x'] < r[['ul_x', 'lr_x']].max()) &
    (df_train['y'] > r[['ul_y', 'lr_y']].min()) &
    (df_train['y'] < r[['ul_y', 'lr_y']].max())]) for i, r in df_sets.iterrows()]
    n_samples = [int(len(df) * .63) for i, df in obs]
    df_sets['n_samples'] = n_samples
    
    # Get support bootstrapped samples and set IDs
    train_dict = {}
    oob_dict = {}
    keep_sets = []
    for i, df in obs:
        if df_sets.ix[i, 'n_samples'] < min_obs:
            continue
        inds = random.sample(df.index, int(len(df) * .63))
        train_dict[i] = inds
        keep_sets.append(i)
        # Get the out of bag samples
        oob_inds = df.index[~df.index.isin(inds)]
        oob_dict[i] = oob_inds
        
    # Drop any sets that don't contain enough training observations
    total_sets = len(df_sets)
    df_drop = df_sets[~df_sets.index.isin(keep_sets)] # Sets not in keep_sets
    n_dropped = len(df_drop)
    print '%s of %s sets dropped because they contained too few observations\n' % (n_dropped, total_sets)
    df_sets = df_sets.ix[keep_sets]
    
    return train_dict, df_sets, oob_dict, df_drop, df_test
    

def coords_to_shp(df, prj_shp, out_shp):
    '''
    Write a shapefile of rectangles from coordinates in df. Each row in df
    represents a unique set of coordinates of a rectangle feature.
    '''
    # Get spatial reference
    ds_prj = ogr.Open(prj_shp)
    lyr_prj = ds_prj.GetLayer()
    srs = lyr_prj.GetSpatialRef()
    ds_prj.Destroy()
    
    # Create output datasource
    driver = ogr.GetDriverByName('ESRI Shapefile')
    #out_shp = os.path.join(out_dir, 'gsrd.shp')
    if os.path.exists(out_shp):
        driver.DeleteDataSource(out_shp)
    out_ds = driver.CreateDataSource(out_shp)
    out_lyr = out_ds.CreateLayer(os.path.basename(out_shp).replace('.shp', ''),
                                 srs,
                                 geom_type=ogr.wkbPolygon)
    
    # Add coord fields
    cols = df.columns
    for c in cols:
        out_lyr.CreateField(ogr.FieldDefn(c, ogr.OFTInteger))
    out_lyr.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    lyr_def = out_lyr.GetLayerDefn()
    
    # Create geometry and add to layer for each row (i.e., feature) in df
    coord_cols = ['ul_x', 'ul_y', 'lr_x', 'lr_y']
    for i, row in df.iterrows():
        
        # Set fields
        feat = ogr.Feature(lyr_def)
        for c in cols:
            feat.SetField(c, row[c])
        feat.SetField('id', i)
        
        # Set geometry
        ul_x, ul_y, lr_x, lr_y = row[coord_cols]
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(ul_x, ul_y) #ul vertex
        ring.AddPoint(lr_x, ul_y) #ur vertex
        ring.AddPoint(lr_x, lr_y) #lr vertex
        ring.AddPoint(ul_x, lr_y) #ll vertex
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        feat.SetGeometry(poly)
        
        # Add the feature to the layer
        out_lyr.CreateFeature(feat)
        feat.Destroy()
    
    out_ds.Destroy()
    
    return out_shp


def get_coords_from_geometry(geom, multipart=False):
    '''Return a list of xy pairs of vertices from a geometry'''
    
    coords = []
    # For each geometry ref, get coordinates from vertices
    for j in xrange(geom.GetGeometryCount()):
        this_g = geom.GetGeometryRef(j) #If geom is a Multipolgyon
        wkt = this_g.ExportToWkt()
        pts_list = wkt.replace('POLYGON','').replace('LINEARRING','').replace('(','').replace(')','').strip().split(',')
        x = [float(p.split(' ')[0]) for p in pts_list]
        y = [float(p.split(' ')[1]) for p in pts_list]
        pts = zip(x,y)
        if multipart == 'split': # Merge multipart features to single part
            coords.extend(pts)
        else:
            coords.append(pts)

    return coords


def get_coords(shp, multipart=None):
    ''' 
    Return a list of lists of the projected coordinates of vertices of shp. 
    Each list represents a feature in the dataset.
    '''    
    ds = ogr.Open(shp)
    if ds == None:
        print 'Shapefile does not exist or is not valid:\n%s' % shp
        return None
    lyr = ds.GetLayer()
    
    # For each feature, get coords
    coords = []
    for i in xrange(lyr.GetFeatureCount()): 
        this_lyr = lyr.GetFeature(i)
        geom = this_lyr.GetGeometryRef()
        these_coords = get_coords_from_geometry(geom)
        if multipart == 'split':
            coords.extend(these_coords)
        else:
            coords.append(these_coords)
    
    # Need the bounds for calculating relative xy in image coords
    bounds = lyr.GetExtent()
    ds.Destroy()

    return coords, bounds


def plot_sets_on_shp(ds_coords, max_size, df_sets, support_size, out_dir=None, fill='0.5', line_color='0.3', line_width=1.0, pad=20, label_sets=False):
    '''
    Plot each rectangle in df_sets overtop of an image of the vector 
    shapes defined by ds_coords. 
    '''
    # Get extreme x's and y's of all sets and compute ratio to rows and cols
    '''minmin_x, maxmax_x = df_sets['ul_x'].min(), df_sets['lr_x'].max()
    minmin_y, maxmax_y = df_sets['lr_y'].min(), df_sets['ul_y'].max()#'''
    maxmax_x = max([max([xy[0] for xy in feature]) for feature in ds_coords]) + support_size[1] * .75
    minmin_x = min([min([xy[0] for xy in feature]) for feature in ds_coords]) - support_size[1] * .75
    maxmax_y = max([max([xy[1] for xy in feature]) for feature in ds_coords]) + support_size[0] * .75
    minmin_y = min([min([xy[1] for xy in feature]) for feature in ds_coords]) - support_size[0] * .75
    delta_x = maxmax_x - minmin_x
    delta_y = maxmax_y - minmin_y
    # Figure out which dimension is larger, make it max_size, and calculate
    #    the proportional size of the other dimension
    if delta_x >= delta_y: 
        cols = max_size
        rows = int(max_size * delta_y/delta_x)
    else:
        cols = int(max_size * delta_x/delta_y)
        rows = max_size
    
    # Create the plot
    fig = plt.figure(figsize=(cols/72.0, rows/72.0), dpi=72)#, tight_layout={'pad':pad})
    sub = fig.add_subplot(1, 1, 1, axisbg='w', frame_on=False)
    sub.axes.get_yaxis().set_visible(False)
    sub.axes.get_xaxis().set_visible(False) 
    sub.set_xlim([minmin_x, maxmax_x])
    sub.set_ylim([minmin_y, maxmax_y])
    
    # Make a list of patches where each feature is a separate patch 
    patches = []
    ds_coords = [c for c in ds_coords if c != None]
    for feature in ds_coords:
        img_coords = [(pt[0], pt[1]) for pt in feature]
        ar_coords = np.array(img_coords)
        poly = matplotlib.patches.Polygon(ar_coords)
        patches.append(poly)
        #sub.add_patch(poly, facecolor=fill, lw=line_width, edgecolor=line_color)
    
    # Make the patch collection and add it to the plot
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, color=fill, lw=line_width, edgecolor=line_color)
    sub.add_collection(p)
    
    # Plot the support sets
    for ind, r in df_sets.iterrows():
        #print ind
        ll_x = r['ul_x']# - x_min) * x_scale
        ll_y = r['lr_y']# - y_min) * y_scale
        w = support_size[1]# * x_scale
        h = support_size[0]# * y_scale
        #if ll_y == minmin_y:
        sub.add_patch(plt.Rectangle((ll_x, ll_y), w, h, facecolor='b', ec='none', lw='0.5', alpha=.1, label=ind))
        if label_sets:
            plt.text(r.ctr_x, r.ctr_y, str(ind), ha='center')

    if out_dir:
        plt.savefig(os.path.join(out_dir, 'support_sets.png'))
    else:
        plt.show()


def tx_from_shp(shp, x_res, y_res):
    ''' Return a GeoTransform of a shapefile if it were rasterized '''
    
    ds = ogr.Open(shp)
    lyr = ds.GetLayer()
    min_x, max_x, min_y, max_y = lyr.GetExtent()
    
    if x_res > 0:
        ul_x = min_x
    else:
        ul_x = max_x
    
    if y_res > 0:
        ul_y = min_y
    else:
        ul_y = max_y
    
    tx = ul_x, x_res, 0, ul_y, 0, y_res
    
    return tx, (min_x, max_x, min_y, max_y)
    
    
def get_gsrd(mosaic_path, cell_size, support_size, n_sets, df_train, min_obs, target_col, predict_cols, out_txt=None, extent_shp=None, pct_train=None, resolution=30, snap_coord=None):
    
    print 'Generating GSRD grid...%s\n' % time.ctime(time.time())
    if mosaic_path.endswith('.shp'):
        x_res = resolution
        y_res = -resolution
        tx, extent = tx_from_shp(mosaic_path, x_res, y_res)
        min_x, max_x, min_y, max_y = extent
    else:
        min_x, min_y, max_x, max_y, x_res, y_res, tx = get_raster_bounds(mosaic_path)
    min_x, min_y, max_x, max_y = [int(c) for c in [min_x, min_y, max_x, max_y]]
    
    # Make the randomly placed grid
    cells = generate_gsrd_grid(cell_size, min_x, min_y, max_x, max_y, x_res, y_res)
    df_grid = pd.DataFrame(cells, columns=['ul_x', 'ul_y', 'lr_x', 'lr_y'])
    out_dir = os.path.dirname(out_txt)
    grid_shp = os.path.join(out_dir, 'gsrd_grid.shp')                         
    coords_to_shp(df_grid, extent_shp, grid_shp)
    print 'Shapefile of support of grid cells written to:\n%s\n' % grid_shp
    
    # Get support sets
    y_size = int(support_size[0]/x_res) * x_res
    x_size = int(support_size[1]/y_res) * y_res
    support_sets = [sample_gsrd_cell(n_sets, bounds, x_size, y_size,
                                     x_res, y_res, tx, snap_coord) for bounds in cells]
    df_sets = pd.concat(support_sets)
    df_sets.index = xrange(len(df_sets))# Original index is range(0,n_sets)
    
    print 'Sampling observations within support sets...'
    t1 = time.time()
    train_dict, df_sets, oob_dict, df_drop, df_test = get_obs_within_sets(df_train, df_sets, min_obs, pct_train)
        
    set_shp = os.path.join(out_dir, 'gsrd_sets.shp')
    coords_to_shp(df_sets, extent_shp, set_shp)
    coords_to_shp(df_drop, extent_shp, set_shp.replace('_sets.shp', '_dropped.shp'))
    print 'Shapefile of support sets written to:\n%s' % set_shp
    print 'Time for sampling: %.1f minutes\n' % ((time.time() - t1)/60)
    
    # Plot the support sets
    '''if extent_shp:
        print 'Plotting support sets...%s\n' % time.ctime(time.time())
        if out_txt: out_dir = os.path.join(os.path.dirname(out_txt))
        else: out_dir = None         
        coords, extent = get_coords(extent_shp)
        set_inds = df_train.set_id.unique()
        plot_sets_on_shp(coords, 900, df_sets.ix[set_inds], support_size, out_dir)'''
    
    # Write train and test dfs to text files
    print 'Writing sample dictionaries to disk...'
    t1 = time.time()
    train_path = os.path.join(out_txt.replace('.txt', '_train_dict.pkl'))
    with open(train_path, 'wb') as f:
        pickle.dump(train_dict, f, protocol=-1)
    
    oob_path = os.path.join(out_txt.replace('.txt', '_oob_dict.pkl'))
    with open(oob_path, 'wb') as f:
        pickle.dump(oob_dict, f, protocol=-1)
    
    if len(df_test) > 0:
        df_test.to_csv(out_txt.replace('.txt', '_test.txt'), sep='\t')
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    print 'Train and test dicts written to:\n', os.path.dirname(out_txt), '\n'
    return train_dict, df_sets, oob_dict


def find_file(basepath, search_str, tsa_str=None, path_filter=None):
    '''
    Return the full path within the directory tree /basepath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    if not os.path.exists(basepath):
        sys.exit('basepath does not exist: \n%s' % basepath)
     
    if tsa_str: 
        bp = os.path.join(basepath, tsa_str)
        if not os.path.isdir(bp):
            exists = False
            # Try adding leading zeros and check for an existing direcoty
            for i in range(10):
                tsa_str = '0' + tsa_str
                bp = os.path.join(basepath, tsa_str)
                if os.path.isdir(bp):
                    exists = True
                    break
            if not exists:
                raise IOError(('Cannot find tile directory with basepath %s and'+\
                'tile_str %s or up to 10 leading zeros') % (basepath, tsa_str))

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
    
    #import pdb; pdb.set_trace()
    if len(paths) < 1:
        #pdb.set_trace()
        sys.exit(('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str))
    
    return paths[0]


def fit_tree_classifier(x_train, y_train, max_features=None):
    ''' Train a decision tree classifier '''
    if not max_features: max_features=None
    dt = tree.DecisionTreeClassifier(max_features=max_features)
    dt.fit(x_train, y_train)
    
    return dt
    
    
def fit_tree_regressor(x_train, y_train, max_features=None):
    ''' Train a decision tree regressor '''
    if not max_features: max_features=None
    dt = tree.DecisionTreeRegressor(max_features=max_features)
    dt.fit(x_train, y_train)
    
    return dt


def fit_bdt_tree_regressor(x_samples, y_samples, pct_bagged=.63):
    ''' 
    Train a decision tree regressor for use in a bagged decision tree.
    The use case is a little different for BDTs than STEM so it's easier just
    to have a separate function.
    '''
    inds = random.sample(x_samples.index, int(len(x_samples) * pct_bagged))
    x_train = x_samples.ix[inds]
    y_train = y_samples.ix[inds]
    dt = tree.DecisionTreeRegressor()
    dt.fit(x_train, y_train)
    
    x_oob = x_samples.ix[~x_samples.index.isin(inds)]
    y_oob = y_samples.ix[x_oob.index]
    oob_rate = calc_oob_rate(dt, y_oob, x_oob)
    
    return dt, oob_rate
    
    
def write_decisiontree(dt, filename):
    ''' 
    Pickle a decision tree and write it to filename. Return the filename.
    '''
    with open(filename, 'w+') as f:
        pickle.dump(dt, f, protocol=-1) #protocol kw might hreduce load time
        
    return filename
    

def write_model(out_dir, df_sets):
    '''
    Write STEM decision trees and dataframe of locations to disk
    '''
    this_dir = os.path.join(out_dir, 'decisiontree_models')
    if not os.path.exists(this_dir):
        os.mkdir(this_dir)
    
    stamp = os.path.basename(out_dir)
    dt_bn = stamp + '_decisiontree_%s'
    
    dt_file = os.path.join(this_dir, dt_bn)
    df_sets['dt_file'] = [write_decisiontree(row.dt_model, dt_file % set_id)\
                         for set_id, row in df_sets.iterrows()]
    
    set_txt = os.path.join(this_dir, stamp + '_support_sets.txt')
    df_sets['set_id'] = df_sets.index
    df_sets.drop('dt_model', axis=1).to_csv(set_txt, sep='\t', index=False)
    
    print 'Support set dataframe and decision trees written to:\n', this_dir
    
    return df_sets, set_txt


def calc_oob_rate(dt, oob_samples, oob_predictors, model_type='classifier'):
    '''
    Return the Out of Bag accuracy rate for the decision tree, dt, and the
    Out of Bag samples, oob_samples
    '''
    
    oob_prediction = dt.predict(oob_predictors)
    
    # If the model is a regressor, calculate the R squared
    if model_type.lower() == 'regressor':
        oob_rate = metrics.r2_score(oob_samples, oob_prediction)
    
    # Otherwise, get the percent correct
    else:
        n_correct = len(oob_samples[oob_samples == oob_prediction])
        n_samples = len(oob_samples)
        oob_rate = int(round(float(n_correct)/n_samples * 100, 0))
    
    return oob_rate


def get_oob_rates(df_sets, df_train, oob_dict, target_col, predict_cols, min_oob=0):
    """
    Calculate the out-of-bag accuracy for each support set in `df_sets`
    
    Parameters
    ----------
    df_sets : pandas DataFrame
        DataFrame with support set info
    df_train : DataFrame of training samples
    oob_dict : dict
        Dictionary of sample indices indicating which samples are OOB for each support set in 
    target_col : string of the target column name
    predict_cols : list-like of string prediction column names
    min_oob : int minimum OOB score. All sets with OOB score < min_oob will be dropped. Default = 0.
    
    Returns
    -------
    df_sets : pandas DataFrame
        Copy of `df_sets` with new columns `oob_rate` appended
    low_oob : pandas DataFrame
        View of df_sets where OOB rate < `min_oob`

    """
    
    for i, row in df_sets.iterrows():
        dt = row.dt_model
        this_oob = df_train.ix[oob_dict[i]]
        oob_samples = this_oob[target_col]
        oob_predictors = this_oob[predict_cols]
        oob_rate = calc_oob_rate(dt, oob_samples, oob_predictors)
        df_sets.ix[i, 'oob_rate'] = oob_rate
    low_oob = df_sets[df_sets.oob_rate < min_oob]
    
    return df_sets, low_oob


def oob_map(ysize, xsize, nodata, nodata_mask, n_tiles, tx, support_size, oob_dict, df_sets, df_train, val_col, var_cols, out_dir, file_stamp, prj, driver):
    
    t0 = time.time()
    
    ul_x, x_res, _, ul_y, _, y_res = tx
    
    df_tiles, df_tiles_rc, tile_size = get_tiles(n_tiles, xsize, ysize, tx)
    total_tiles = len(df_tiles)
    #coords_to_shp(df_tiles, '/vol/v2/stem/conus_testing/models/landcover_20170420_1347/gsrd_sets.shp', '/vol/v2/stem/conus_testing/models/landcover_20170503_1820/ag_tiles.shp')
    # Find the tiles that have only nodata values
    '''t1 = time.time()
    print '\nFinding empty tiles...'
    empty_tiles = find_empty_tiles(df_tiles, nodata_mask, tx)
    # Select only tiles that are not empty
    df_tiles = df_tiles.select(lambda x: x not in empty_tiles)
    total_tiles = len(df_tiles)
    print '%s empty tiles found of %s total tiles\n' %\
    (len(empty_tiles), total_tiles)'''
    
    # Make directory for storing tiles 
    tile_dir = os.path.join(out_dir, 'oob_tiles_temp')
    if not os.path.exists(tile_dir):
        os.mkdir(tile_dir)
    
    if 'oob_rate' not in df_sets.columns:
        for i, row in df_sets.iterrows():
            with open(row.dt_file) as f:
                dt = pickle.load(f)
            this_oob = df_train.ix[oob_dict[i]]
            oob_samples = this_oob[val_col]
            oob_predictors = this_oob[var_cols]
            df_sets.ix[i, 'oob_rate'] = calc_oob_rate(dt, oob_samples, oob_predictors)#'''
    empty_tiles = []
    for i, (tile_id, t_row) in enumerate(df_tiles[4:].iterrows()):
        t1 = time.time()
        print 'Aggregating for %s of %s tiles' % (i + 1, total_tiles)
        print 'Tile index: ', tile_id
        
        rc = df_tiles_rc.ix[tile_id]
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[tile_id]
        
        # Get a mask for this tile
        if type(nodata_mask) == np.ndarray:
            tile_mask = nodata_mask[ul_r : lr_r, ul_c : lr_c]
        else: # Otherwise, it's a path to a shapefile
            tile_mask, _ = mosaic_by_tsa.kernel_from_shp(nodata_mask, t_row[['ul_x', 'ul_y', 'lr_x', 'lr_y']], tx, -9999)
            tile_mask = tile_mask != -9999
        # If it's empty, skip and add it to the list
        if not tile_mask.any():
            empty_tiles.append(tile_id)
            print 'Tile empty. Skipping...\n'
            continue

        # Calculate the size of this tile in case it's at the edge where the
        #   tile size will be slightly different
        this_size_xy = abs(t_row.lr_y - t_row.ul_y), abs(t_row.lr_x - t_row.ul_x)
        df_these_sets = get_overlapping_sets(df_sets, t_row, this_size_xy, support_size)
        
        this_size = rc.lr_r - rc.ul_r, rc.lr_c - rc.ul_c
        n_sets = len(df_these_sets)
        tile_ul = t_row[['ul_x','ul_y']]
        
        print n_sets, ' Overlapping sets'
        t2 = time.time()
        oob_rates = []
        oob_bands = []
        for s_ind, s_row in df_these_sets.iterrows():
            #set_coords = s_row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
            # Fill a band for this array
            offset = calc_offset(tile_ul, (s_row.ul_x, s_row.ul_y), tx)
            tile_inds, a_inds = mosaic_by_tsa.get_offset_array_indices(tile_size, (support_size[0]/30, support_size[1]/30), offset)
            nrows = a_inds[1] - a_inds[0]
            ncols = a_inds[3] - a_inds[2]
            try:
                ar_oob_rate = np.full((nrows, ncols), s_row.oob_rate, dtype=np.int16)
            except:
                import pdb; pdb.set_trace()
            
            oob_band = np.full(this_size, nodata)
            oob_band[tile_inds[0]:tile_inds[1], tile_inds[2]:tile_inds[3]] = ar_oob_rate
            oob_band = oob_band.astype(np.float16)
            oob_band[~tile_mask | (oob_band==nodata)] = np.nan
            oob_bands.append(oob_band)
            oob_rates.append(s_row.oob_rate)
         
        print 'Average OOB: ', int(np.mean(oob_rates))
        print 'Filling tiles: %.1f seconds' % ((time.time() - t2))
        ar_tile = np.dstack(oob_bands)
        del oob_bands
        t2 = time.time()
        this_oob = np.nanmean(ar_tile, axis=2).astype(np.int16)
        this_cnt = np.sum(~np.isnan(ar_tile), axis=2).astype(np.int16)
        
        nans = np.isnan(this_oob)
        this_oob[nans] = nodata
        this_cnt[nans] = nodata
        path_template = os.path.join(tile_dir, 'tile_%s_%s.tif')
        this_tx = t_row.ul_x, tx[1], tx[2], t_row.ul_y, tx[4], tx[5]
        mosaic_by_tsa.array_to_raster(this_oob, this_tx, prj, driver, path_template % (tile_id, 'oob'), nodata=nodata, silent=True)
        mosaic_by_tsa.array_to_raster(this_cnt, this_tx, prj, driver, path_template % (tile_id, 'count'), nodata=nodata, silent=True)
        print 'Aggregating: %.1f minutes' % ((time.time() - t2)/60)
        print 'Total time for this tile: %.1f minutes\n' % ((time.time() - t1)/60)
        
    # For count and oob, stitch tiles back together
    avg_dict = {}
    for stat in ['oob', 'count']: 
        ar = np.full((ysize, xsize), nodata, dtype=np.uint8)
        for tile_id, tile_coords in df_tiles.iterrows():
            if tile_id in empty_tiles:
                continue
            rc = df_tiles_rc.ix[tile_id]
            tile_rows = rc.lr_r - rc.ul_r
            tile_cols = rc.lr_c - rc.ul_c
            tile_file = os.path.join(tile_dir, 'tile_%s_%s.tif' % (tile_id, stat))
            if not os.path.isfile(tile_file):
                continue
            ds = gdal.Open(tile_file)
            ar_tile = ds.ReadAsArray()
            t_ulx = tile_coords[['ul_x', 'ul_y']]
            row_off, col_off = calc_offset((ul_x, ul_y), t_ulx, tx)
            try:
                ar[row_off : row_off + tile_rows, col_off : col_off + tile_cols] = ar_tile
            except:
                import pdb; pdb.set_trace()
        out_path = os.path.join(out_dir, '%s_%s.tif' % (file_stamp, stat))
        gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(ar.dtype)
        mosaic_by_tsa.array_to_raster(ar, tx, prj, driver, out_path, gdal_dtype, nodata=nodata)
        avg_dict[stat] = ar[ar != nodata].mean()
    print '\nTotal aggregation run time: %.1f hours' % ((time.time() - t0)/3600)
    
    # Delete the temporary tile directory
    shutil.rmtree(tile_dir)
    
    return avg_dict, df_sets#"""


def get_predict_array(args):

    ar = mosaic_by_tsa.get_mosaic(*args[1:])
    return args[0], ar.ravel()
    
    
def get_predictors(df_var, mosaic_tx, tile_strs, tile_ar, coords, nodata_mask, out_nodata, set_id, constant_vars=None):
    '''
    Return an array of flattened predictor arrays where each predictor is a 
    separate column
    '''
    t0 = time.time()
    predictors = []
    for ind, var in enumerate(df_var.index):
        #this_tile_ar = np.copy(tile_ar)
        data_band, search_str, basepath, by_tile, path_filter = df_var.ix[var]
        if constant_vars: search_str = search_str.format(constant_vars['YEAR'])
        
        if by_tile:
            files = [find_file(basepath, search_str, tile, path_filter) for tile in tile_strs]
            ar_var = mosaic_by_tsa.get_mosaic(mosaic_tx, tile_strs, tile_ar, coords, data_band, set_id, files)
            
        else:
            try:
                this_file = find_file(basepath, search_str, path_filter=path_filter)
                tx_, ar_var, roff_, coff_ = mosaic_by_tsa.get_array(this_file, data_band, coords)
            except:
                exc_type, exc_msg, _ = sys.exc_info()
                exc_type_str = str(exc_type).split('.')[-1].replace("'", '').replace('>', '')
                import pdb; pdb.set_trace()
        ar_var[nodata_mask] = out_nodata
        predictors.append(ar_var.ravel())

    #df_predict = pd.DataFrame(predictors)
    if constant_vars:
        size = predictors[0].size
        for const in sorted(constant_vars.keys()):
            val = constant_vars[const]
            predictors.append(np.full(size, val, dtype=np.int16))
   
    ar = np.vstack(predictors).T # each pixel is a row with each predictor in a col
    del predictors
          
    return ar
    


def predict_set(set_id, df_var, mosaic_ds, coords, mosaic_tx, xsize, ysize, dt, nodata, dtype=np.int16, constant_vars=None):
    '''
    Return a predicted array for set with id==set_id
    '''
    # Get an array of tile_ids within the bounds of coords
    if type(mosaic_ds) == ogr.Layer:
        tile_ar, tile_off = mosaic_by_tsa.kernel_from_shp(mosaic_ds, coords, mosaic_tx, nodata)
    else:
        tile_ar, tile_off = mosaic_by_tsa.extract_kernel(mosaic_ds, 1, coords, mosaic_tx,
                                                xsize, ysize, nodata=nodata)
    tile_mask = tile_ar == 0
    tile_ar[tile_mask] = nodata
    
    # Get the ids of TSAs this kernel covers
    tile_ids = np.unique(tile_ar)
    tile_strs = [str(t) for t in tile_ids if t != nodata]
    array_shape = tile_ar.shape

    # Get an array of predictors where each column is a flattened 2D array of a
    #   single predictor variable
    temp_nodata = -9999
    ar_predict = get_predictors(df_var, mosaic_tx, tile_strs, tile_ar, coords, tile_mask, temp_nodata, set_id, constant_vars)
    del tile_ar #Release resources from the tsa array
    
    t0 = time.time()
    nodata_mask = ~ np.any(ar_predict==temp_nodata, axis=1)
    predictions = dt.predict(ar_predict[nodata_mask]).astype(dtype)
    ar_prediction = np.full(ar_predict.shape[0], nodata, dtype=dtype)
    ar_prediction[nodata_mask] = predictions
    
    #print 'Time for predicting: %.1f minutes' % ((time.time() - t0)/60)
    
    return ar_prediction.reshape(array_shape)


def par_predict(args):
    '''Helper function to parallelize predicting with multiple decision trees'''
    
    t0 = time.time()

    coords, mosaic_type, mosaic_path, mosaic_tx, prj, nodata, set_count, total_sets, set_id, df_var, xsize, ysize, dt_file, nodata, dtype, constant_vars, predict_dir = args
    print '\nPredicting for set %s of %s' % (set_count + 1, total_sets)
    # Save rasters of tsa arrays ahead of time to avoid needing to pickle or fork mosaic
    if mosaic_type == 'vector':
        mosaic_ds = ogr.Open(mosaic_path)
        mosaic_lyr = mosaic_ds.GetLayer()
        tile_ar, tile_off = mosaic_by_tsa.kernel_from_shp(mosaic_lyr, coords, mosaic_tx, nodata=0)
    else:
        mosaic_ds = gdal.Open(mosaic_path)
        tile_ar, tile_off = mosaic_by_tsa.extract_kernel(mosaic_ds, 1, coords,
                                                mosaic_tx, xsize, ysize,
                                                nodata=nodata)
    mosaic_ds = None
    tx_out = coords.ul_x, mosaic_tx[1], mosaic_tx[2], coords.ul_y, mosaic_tx[4], mosaic_tx[5]
    
    with open(dt_file, 'rb') as f: 
        dt_model = pickle.load(f)
    try:
        tile_mask = tile_ar == 0
        tile_ar[tile_mask] = nodata
        
        # Get the ids of tiles this kernel covers
        tile_ids = np.unique(tile_ar)
        tile_strs = [str(tile) for tile in tile_ids if tile != nodata]
        array_shape = tile_ar.shape
    
        # Get an array of predictors where each column is a flattened 2D array of a
        #   single predictor variable
        temp_nodata = -9999
        ar_predict = get_predictors(df_var, mosaic_tx, tile_strs, tile_ar, coords, tile_mask, temp_nodata, set_id, constant_vars)
        del tile_ar
        
        # Predict only for pixels where no predictors are nodata
        nodata_mask = ~ np.any(ar_predict==temp_nodata, axis=1)  
        predictions = dt_model.predict(ar_predict[nodata_mask]).astype(dtype)
        ar_prediction = np.full(ar_predict.shape[0], nodata, dtype=dtype)
        ar_prediction[nodata_mask] = predictions
        
        # Save the predicted raster
        driver = gdal.GetDriverByName('gtiff')
        ar_prediction = ar_prediction.reshape(array_shape)
        out_path = os.path.join(predict_dir, 'prediction_%s.tif' % set_id)
        np_dtype = get_min_numpy_dtype(ar_prediction)
        gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(np_dtype)
        mosaic_by_tsa.array_to_raster(ar_prediction, tx_out, prj, driver, out_path, gdal_dtype, nodata=nodata)
        
    except:
        print 'problem with set_count: ', set_count, 'set_id: ', set_id
        errorLog = os.path.dirname(predict_dir)+'/predication_errors.txt'       
        if not os.path.isfile(errorLog):                 
            with open(errorLog, 'w') as el:       
                el.write('set_count: '+str(set_count)+'\n')
                el.write('set_id: '+str(set_id)+'\n')
                el.write(traceback.format_exc()+'\n')
        else:
            with open(errorLog, 'a') as el:       
                el.write('set_count: '+str(set_count)+'\n')
                el.write('set_id: '+str(set_id)+'\n')
                el.write(traceback.format_exc()+'\n')
        return
        #sys.exit(traceback.print_exception(*sys.exc_info()))
    
    print 'Total time for set %s of %s: %.1f minutes' % (set_count + 1, total_sets, (time.time() - t0)/60)
        

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


def get_overlapping_sets(df_sets, tile_bounds, tile_size, support_size, resolution=30):
    '''
    Return a dataframe of support sets that overlap the tile defined by
    tile bounds
    '''
    set_ysize, set_xsize = support_size
    tile_rows, tile_cols = tile_size
    
    # Calculate the max distance that the centers of the tile and the set
    #   could be from one another
    max_x_dist = (set_xsize/2 + tile_cols/2)/resolution
    max_y_dist = (set_ysize/2 + tile_rows/2)/resolution
    
    # Get sets where both the x and y center coords are within the respective
    #   max distances of one another
    x_dif = (df_sets.ctr_x - tile_bounds.ctr_x).abs().astype(int)/resolution
    y_dif = (df_sets.ctr_y - tile_bounds.ctr_y).abs().astype(int)/resolution 
    overlap = df_sets[(x_dif < max_x_dist) & (y_dif < max_y_dist)]
    #import pdb; pdb.set_trace()
    
    return overlap


def calc_offset(ul_xy1, ul_xy2, tx):
    '''
    Return the row and col offset of a data array from a tsa_array
    '''
    x1, y1 = ul_xy1
    x2, y2 = ul_xy2
    row_off = int((y2 - y1)/tx[5])
    col_off = int((x2 - x1)/tx[1])
    
    #return pd.Series((row_off, col_off))
    return row_off, col_off

def shp_to_df(shp):
    ''' Return a pandas dataframe of a shapefile attribute table'''    

    ds = ogr.Open(shp)
    lyr = ds.GetLayer()
    lyr_def = lyr.GetLayerDefn()
    
    fields = [lyr_def.GetFieldDefn(i).GetName() for i in range(lyr_def.GetFieldCount())]
    attributes = []
    for i in xrange(lyr.GetFeatureCount()):
        feat = lyr.GetFeature(i)
        these_vals = {field: feat.GetField(field) for field in fields}
        these_vals['fid'] = feat.GetFID()
        attributes.append(these_vals)
    
    lyr.Destroy()
    ds.Destroy()
    
    df = pd.DataFrame(attributes)
    df.set_index('fid', inplace=True)
    
    return df
    
def load_predictions(p_dir, df_sets, tile_ul, tile_size):

    predictions = {}
    
    for set_id in df_sets.index:
        f = os.path.join(p_dir, 'prediction_%s.tif' % set_id)
        ds = gdal.Open(f)
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        tx = ds.GetGeoTransform()
        offset = calc_offset(tile_ul, (tx[0], tx[3]), tx)
        t_inds, a_inds = mosaic_by_tsa.get_offset_array_indices(tile_size, (ysize, xsize), offset)
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
        ar_pred = ar_pred.astype(np.float16)
        ar_pred[ar_pred == nodata] = np.nan

        ar_tile[tile_inds[0]:tile_inds[1], tile_inds[2]:tile_inds[3]] = ar_pred
        #ar_pred[set_row_u:set_row_d, set_col_l:set_col_r]
    except Exception as e:
        #import pdb; pdb.set_trace()
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
    ar_out = np.full(ar.shape, nodata, dtype=np.int16)
    data_mask = ar != nodata
    ar_out[data_mask] = mp[ar[data_mask]]
    
    return ar_out
  
  
def find_empty_tiles(df, nodata_mask, tx, nodata=None):
    
    empty = []
    
    for i, coords in df.iterrows():
        ul_r, ul_c = calc_offset((tx[0], tx[3]), coords[['ul_x', 'ul_y']], tx)
        lr_r, lr_c = calc_offset((tx[0], tx[3]), coords[['lr_x', 'lr_y']], tx)
        if type(nodata_mask) == ogr.Layer:
            this_mask, _ = mosaic_by_tsa.kernel_from_shp(nodata_mask, coords, tx, nodata)
        else:
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
    
    return np.round(w_mean,0).astype(np.int16)
    

def pct_vote(ar, ar_vote, ar_count):
    
    shape = ar.shape
    ar_eq = ar == ar_vote.repeat(shape[2]).reshape(shape)
    ar_sum = ar_eq.sum(axis=2)
    ar_pct = np.round(ar_sum/ar_count.astype(np.float16) * 100).astype(np.uint8)
    
    return ar_pct


def aggregate_tile(tile_coords, n_tiles, nodata_mask, support_size, agg_stats, prediction_dir, df_sets, nodata, out_dir, file_stamp, tx, prj, ar_tile=None):
    
    ############################################################################################################################
    # jdb added 6/22/2017 
    # for testing purposes, skip tiles that have already been aggregated
    path_template = os.path.join(out_dir, 'tile_{0}_*'.format(tile_coords.name))
    #pdb.set_trace()    
    if len(glob.glob(path_template)) != 0:
      print 'Aggregating for %s of %s tiles is already done - skipping...' % (tile_coords.name + 1, n_tiles)
      return
    
    t0 = time.time()
    print 'Aggregating for %s of %s tiles...' % (tile_coords.name + 1, n_tiles)
    ############################################################################################################################
    
    # Get overlapping sets
    x_res = tx[1]
    y_res = tx[5]
    tile_rows = (tile_coords.lr_y - tile_coords.ul_y) / y_res
    tile_cols = (tile_coords.lr_x - tile_coords.ul_x) / x_res
    tile_size = tile_rows, tile_cols
    df_these_sets = get_overlapping_sets(df_sets, tile_coords, tile_size, support_size)
    n_sets = len(df_these_sets)
    
    # Load overlapping predictions from disk and read them as arrays
    tile_ul = tile_coords[['ul_x','ul_y']]
    predictions = load_predictions(prediction_dir, df_these_sets, tile_ul, tile_size)
    
    print n_sets, ' Overlapping sets'
    pred_bands = []
    importance_bands = []
    for s_ind in df_these_sets.index:
        
        # Fill tile with prediction
        ar_pred, tile_inds = predictions[s_ind]
        pred_band = fill_tile_band(tile_size, ar_pred, tile_inds, nodata)
        pred_bands.append(pred_band)
        
        # Get feature with maximum importance and fill tile with that val
        if 'importance' in agg_stats:
            try:
                ar_import = np.full(ar_pred.shape, df_sets.ix[s_ind, 'max_importance'], dtype=np.uint8)
                import_band = fill_tile_band(tile_size, ar_import, tile_inds, nodata)
                importance_bands.append(import_band)
            except Exception as e:
                print e
                continue#'''    
    ar_tile = np.dstack(pred_bands)
    if 'importance' in agg_stats:
        ar_impr = np.dstack(importance_bands)
    del pred_bands, importance_bands#, ar_import
    #print 'Filling tiles: %.1f seconds' % ((time.time() - t1))
    
    driver = gdal.GetDriverByName('gtiff')
    path_template = os.path.join(out_dir, 'tile_{0}_%s.tif'.format(tile_coords.name))
    nans = np.isnan(ar_tile).all(axis=2)
    if 'mean' in agg_stats:
        ar_mean = np.nanmean(ar_tile, axis=2)
        ar_mean[nans | nodata_mask] = nodata
        out_path = path_template %  'mean'
        mosaic_by_tsa.array_to_raster(ar_mean, tx, prj, driver, out_path, nodata=nodata, silent=True)
        ar_mean = None
        
    if 'stdv' in agg_stats or 'stdv' in agg_stats:
        ar_stdv = np.nanstdv(ar_tile, axis=2)
        ar_stdv[nans | ~nodata_mask] = -9999
        out_path = path_template %  'stdv'
        mosaic_by_tsa.array_to_raster(ar_stdv, tx, prj, driver, out_path, nodata=-9999, silent=True)
        ar_stdv = None

    if 'vote' in agg_stats:
        ar_vote = mode(ar_tile, axis=2)
        ar_vote[nans | ~nodata_mask] = nodata
        out_path = path_template %  'vote'
        mosaic_by_tsa.array_to_raster(ar_vote, tx, prj, driver, out_path, nodata=nodata, silent=True)

    if 'count' in agg_stats:
        ar_count = np.sum(~np.isnan(ar_tile), axis=2)
        ar_count[nans | ~nodata_mask] = nodata
        out_path = path_template %  'count'
        mosaic_by_tsa.array_to_raster(ar_count, tx, prj, driver, out_path, nodata=nodata, silent=True)

    if 'importance' in agg_stats:
        ar_impor = mode(ar_impr, axis=2)
        ar_impor[nans | ~nodata_mask] = nodata
        out_path = path_template %  'importance'
        mosaic_by_tsa.array_to_raster(ar_impor, tx, prj, driver, out_path, nodata=nodata, silent=True)
        ar_impor = None
            
    if 'pct_vote' in agg_stats:
        if not 'vote' in agg_stats:
            ar_vote = np.mode(ar_tile, axis=2)
        if not 'count' in agg_stats:
            ar_count = np.sum(~np.isnan(ar_tile), axis=2)  
        ar_pcvt = pct_vote(ar_tile, ar_vote, ar_count)
        ar_pcvt[nans | ~nodata_mask] = nodata
        out_path = path_template %  'pct_vote'
        mosaic_by_tsa.array_to_raster(ar_pcvt, tx, prj, driver, out_path, nodata=nodata, silent=True)
    
    if 'oob' in agg_stats:
        ar_oob = np.nanmean(ar_tile, axis=2)
        ar_oob[nans | ~nodata_mask] = nodata
        out_path = path_template %  'oob'
        mosaic_by_tsa.array_to_raster(ar_oob, tx, prj, driver, out_path, nodata=nodata, silent=True)
        ar_oob = None
    
    print 'Total time for tile %s: %.1f minutes\n' % (tile_coords.name, ((time.time() - t0)/60)) 


def par_aggregate_tile(args):
    
    aggregate_tile(*args)


def aggregate_predictions(n_tiles, ysize, xsize, nodata, nodata_mask, tx, support_size, agg_stats, prediction_dir, df_sets, out_dir, file_stamp, prj, driver, n_jobs=None):
    
    #mosaic_ds = gdal.Open(mosaic_path)
    df_tiles, df_tiles_rc, tile_size = get_tiles(n_tiles, xsize, ysize, tx)
    #coords_to_shp(df_tiles, '/vol/v2/stem/conus_testing/models/landcover_20170419_1443/gsrd_sets.shp', '/vol/v2/stem/conus_testing/models/landcover_20170419_1443/ag_tiles.shp')
    n_tiles = len(df_tiles)
    df_tiles['tile'] = df_tiles.index
    
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
    pct_importance = importance / importance.sum()
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    
    tile_dir = os.path.join(out_dir, 'tiles')
    if not os.path.exists(tile_dir):
        os.mkdir(tile_dir)
    # Loop through each tile
    args = []
    empty_tiles = []
    for ind, tile_coords in df_tiles.iterrows():
        tx_tile = tile_coords.ul_x, x_res, x_rot, tile_coords.ul_y, y_rot, y_res
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[ind]
        # if running in parallel, add to list of args
        # Get a mask for this tile
        if type(nodata_mask) == np.ndarray:
            tile_mask = nodata_mask[ul_r : lr_r, ul_c : lr_c]
        else: # Otherwise, it's a path to a shapefile
            tile_mask, _ = mosaic_by_tsa.kernel_from_shp(nodata_mask, tile_coords[['ul_x', 'ul_y', 'lr_x', 'lr_y']], tx_tile, -9999)
            #mosaic_by_tsa.array_to_raster(tile_mask, tx_tile, prj, driver, os.path.join(tile_dir, 'delete_%s.tif' % ind), nodata=-9999)
            tile_mask = tile_mask != -9999
           
        # If it's empty, skip and add it to the list
        if not tile_mask.any():
            empty_tiles.append(ind)
            print 'Skipping tile %s of %s because it contains only nodata values...\n' % (ind + 1, n_tiles)
            continue
        if n_jobs:
            args.append([tile_coords, n_tiles, tile_mask, support_size, agg_stats, prediction_dir,
                         df_sets, nodata, tile_dir, file_stamp, tx_tile, prj])
        # Otherwise, just aggregate for this tile
        else:
            try:
                aggregate_tile(tile_coords, n_tiles, tile_mask, support_size, agg_stats, prediction_dir,
                           df_sets, nodata, tile_dir, file_stamp, tx_tile, prj)
            except:
                print 'Problem with tile ', ind
                print traceback.print_exception(*sys.exc_info())
    
    # agregate in parallel if n_jobs is given
    if n_jobs:
        p = Pool(n_jobs)
        p.map(par_aggregate_tile, args, 1)#'''
    
    # For each stat given, load each tile, add it to an array, and save the array
    #xsize = mosaic_ds.
    driver = gdal.GetDriverByName('gtiff')
    for stat in agg_stats:
        if stat == 'stdv':
            this_nodata = -9999
            ar = np.full((ysize, xsize), this_nodata, dtype=np.int16) 
        else:
            this_nodata = nodata
            ar = np.full((ysize, xsize), this_nodata, dtype=np.uint8)
        
        for tile_id, tile_coords in df_tiles.iterrows():
            if tile_id in empty_tiles:
                continue
            tile_rows = (tile_coords.lr_y - tile_coords.ul_y) / y_res
            tile_cols = (tile_coords.lr_x - tile_coords.ul_x) / x_res
            tile_file = os.path.join(tile_dir, 'tile_%s_%s.tif' % (tile_id, stat))
            ds = gdal.Open(tile_file)
            ar_tile = ds.ReadAsArray()
            t_ulx = tile_coords[['ul_x', 'ul_y']]
            row_off, col_off = calc_offset((ul_x, ul_y), t_ulx, tx)
            ar[row_off : row_off + tile_rows, col_off : col_off + tile_cols] = ar_tile
        
        out_path = os.path.join(out_dir, '%s_%s.tif' % (file_stamp, stat))
        gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(ar.dtype)
        mosaic_by_tsa.array_to_raster(ar, tx, prj, driver, out_path, gdal_dtype, nodata=this_nodata)
    
    # Clean up the tiles
    shutil.rmtree(tile_dir)
    
    return pct_importance, df_sets


def evaluate_ebird(sample_txt, ar, tx, cell_size, target_col, n_per_cell, n_trials=50, year=None):
    t0 = time.time()
    
    df_test = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    #df_test = pd.concat([df_test[df_test[target_col] == 1].drop_duplicates(subset=['row', 'col']), df_test[df_test[target_col] == 0]])
    #df_test.drop_duplicates(subset=[target_col, 'row', 'col'], inplace=True)

    xsize, ysize = ar.shape
    df_test['predicted'] = ar[df_test.row, df_test.col]
    df_test = df_test[df_test.predicted != 255]
    if year:
        df_test = df_test[df_test.YEAR == year]
    
    # Get bounds for cells in a GSRD
    ul_x, x_res, _, ul_y, _, y_res = tx
    lr_x = xsize * x_res/abs(x_res)
    lr_y = ysize * y_res/abs(y_res)
    min_x = min([ul_x, lr_x])
    max_x = max([ul_x, lr_x])
    min_y = min([ul_y, lr_y])
    max_y = max([ul_y, lr_y])
    cells = generate_gsrd_grid(cell_size, min_x, min_y, max_x, max_y, x_res, y_res)
    
    
    # Get all unique sample locations within each cell
    print 'Getting samples within each cell...'
    t1 = time.time()
    locations = []
    for i, (ul_x, ul_y, lr_x, lr_y) in enumerate(cells):
        # Get all samples within this cell
        df_temp = df_test[
                        (df_test.x > min([ul_x, lr_x])) &
                        (df_test.x < max([ul_x, lr_x])) &
                        (df_test.y > min([ul_y, lr_y])) &
                        (df_test.y > max([ul_y, lr_y]))
                        ]
        # Get all unique location
        unique = [(rc[0],rc[1]) for rc in df_temp[['row','col']].drop_duplicates().values]
        locations.append([i, np.array(unique)])
        #locations.append([i, df_temp.index])
    print '%.1f seconds\n' % (time.time() - t1)
    
    # Run n_trials from random samples
    print 'Calculating accuracy for %s trials...' % n_trials
    t1 = time.time()
    results = []
    used_idx = []
    roc_curves = []
    for i in range(n_trials):
        # Get n_per_cell random samples from each cell
        random_locations = []
        for i, l in locations:
            if len(l) < n_per_cell: 
                continue
            #random_locations.extend(l[random.sample(range(len(l)), n_per_cell)])
            random_locations.extend(random.sample(l, n_per_cell))
        rows, cols = zip(*random_locations)
        #idx = []
        #for row, col in random_locations:
        #    idx.extend(df_test[(df_test.row == row) & (df_test.col == col)].index.tolist())
        df_samples = df_test[df_test.row.isin(rows) & df_test.col.isin(cols)]
        #df_samples = df_test.ix[random_locations]
        used_idx.extend(df_samples.index.tolist())
        t_vals = df_samples[target_col]
        p_vals = df_samples.predicted/100.0
        
        # Calc rmspe, rmspe for positive samples, rmspe for neg. samples, and auc
        try:
            auc = round(metrics.roc_auc_score(t_vals, p_vals), 3)
            this_roc_curve = metrics.roc_curve(t_vals, p_vals)
            roc_curves.append(this_roc_curve)
        except:
            auc = 0
        r2 = metrics.r2_score(t_vals, p_vals)
        ac, ac_s, ac_u, ssd, spod = ev.calc_agree_coef(t_vals, p_vals, t_vals.mean(), p_vals.mean())
        rmse = ev.calc_rmse(t_vals, p_vals)
        false_mask = t_vals == 0
        rmse_n = ev.calc_rmse(t_vals[false_mask], p_vals[false_mask])
        rmse_p = ev.calc_rmse(t_vals[~false_mask], p_vals[~false_mask])
        
        
        results.append({'ac': ac,
                        'ac_s': ac_s,
                        'ac_u': ac_u,
                        'r2': r2,
                        'rmse': rmse,
                        'rmse_n': rmse_n,
                        'rmse_p': rmse_p,
                        'auc': auc
                       })
    print '%.1f seconds\n' % (time.time() - t1)
                       
    df = pd.DataFrame(results)
    df_samples = df_test.ix[np.unique(used_idx)]
    
    return df, df_samples, roc_curves
    

def evaluate_by_lc(df, ar, lc_path, target_col, lc_classes=None, ar_nodata=255):
    
    #df = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    
    ds = gdal.Open(lc_path)
    ar_lc = ds.ReadAsArray()
    ds = None
    
    df['predicted'] = ar[df.row, df.col]/100.0
    df = df[df.predicted != ar_nodata]
    df['lc_class'] = ar_lc[df.row, df.col]
    if not lc_classes:
        lc_classes = df.lc_class.unique()
    
    df_stats = pd.DataFrame(columns=['auc','rmse', 'lc_class'])
    for lc in lc_classes:
        df_lc = df[df.lc_class == lc]
        if len(df_lc) == 0:
            print '\nNo samples found in land cover class %s' % lc
            continue
        t_vals = df_lc[target_col]
        p_vals = df_lc.predicted
        if t_vals.min() == t_vals.max():
            print '\nOnly one class present for class %s. Skipping...\n' % lc
            continue
        n_pos = len(t_vals[t_vals == 1])
        n_neg = len(t_vals[t_vals == 0])
        auc = metrics.roc_auc_score(t_vals, p_vals)
        rmse = ev.calc_rmse(t_vals, p_vals)
        lc_dict = {'lc_class': lc, 'rmse': rmse,'auc': auc, 'n_pos': n_pos, 'n_neg': n_neg}
        df_stats = df_stats.append(pd.DataFrame([lc_dict],
                                                 index=[lc]))
    df_stats.set_index('lc_class', inplace=True)
    return df_stats.sort_index()


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
    
    coords = this_set[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
    mosaic_dir = '/vol/v2/stem/canopy/canopy_20160212_2016/var_mosaics'
    saving_stuff = set_id, mosaic_dir, prj, driver
    ar_predict = predict_set(set_id, df_var, mosaic_ds, coords, mosaic_tx, xsize, ysize, dt_model, saving_stuff)
    return ar_predict
    
    '''out_dir = '/vol/v2/stem/scripts/testing'
    out_path = os.path.join(out_dir, 'predict_rerun_%s.bsq' % set_id)
    
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    tx = this_set.ul_x, x_res, x_rot, this_set.ul_y, y_rot, y_res
    mosaic_by_tsa.array_to_raster(ar_predict, tx, prj, driver, out_path, GDT_Int32)'''

def get_gdal_dtype(type_code):
    
    code_dict = {1: gdal.GDT_Byte,
                 2: gdal.GDT_UInt16,
                 3: gdal.GDT_Int16,
                 4: gdal.GDT_UInt32,
                 5: gdal.GDT_Int32,
                 6: gdal.GDT_Float32,
                 7: gdal.GDT_Float64,
                 8: gdal.GDT_CInt16,
                 9: gdal.GDT_CInt32,
                 10: gdal.GDT_CFloat32,
                 11: gdal.GDT_CFloat64
                 }
    
    return code_dict[type_code]




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
#target_col = 'canopy'
#n_tiles = 10
#out_dir = '/vol/v2/stem/canopy/models/'
#set_id, ar, df_sets, df_train = main(sample_txt, target_col, mosaic_path, cell_size, support_size, sets_per_cell, min_obs, pct_train, target_col, n_tiles, out_dir)

#params = '/vol/v2/stem/param_files/build_stem_params_nomse.txt'
#predictions, df_sets, df_train = main(params)
#stuff = main(params)
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
