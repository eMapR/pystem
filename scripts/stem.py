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
import shutil
import random
import warnings
#from itertools import count as itertoolscount
#from random import sample as randomsample
#from string import ascii_lowercase
from osgeo import gdal, ogr
from gdalconst import *
from sklearn import tree, metrics
from multiprocessing import Pool
from datetime import datetime
import cPickle as pickle
import pandas as pd
import numpy as np

# Import ancillary scripts
import mosaic_by_tsa as mosaic
import evaluation as ev

gdal.UseExceptions()
warnings.filterwarnings('ignore')
# Turn off the annoying setting value on a copy warning
pd.options.mode.chained_assignment = None

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
    
    # Need to know if x increases l or r and if y increases up or down
    #x_res_sign = int(x_res/abs(x_res))
    #y_res_sign = int(y_res/abs(y_res))
    
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
    
    # Calculate the ul of each cell
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


def sample_gsrd_cell(n, cell_bounds, x_size, y_size, x_res, y_res, tx):
    '''
    Return a list of bounding coordinates for n support sets from 
    randomly sampled x,y center coords within bounds
    '''
    ul_x, ul_y, lr_x, lr_y = cell_bounds
    min_x, max_x = min(ul_x, lr_x), max(ul_x, lr_x)
    min_y, max_y = min(ul_y, lr_y), max(ul_y, lr_y)
    
    # Calculate support set centers and snap them to the ultimate raster grid
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

    #numeric_bounds = zip(min_x_ls, min_y_ls, max_x_ls, max_y_ls)
    these_bounds = zip(ul_x_ls, ul_y_ls, lr_x_ls, lr_y_ls, x_centers, y_centers)
    #df = pd.DataFrame(numeric_bounds, columns=['min_x', 'min_y', 'max_x', 'max_y'])
    columns = ['ul_x', 'ul_y', 'lr_x', 'lr_y', 'ctr_x', 'ctr_y']
    df = pd.DataFrame(these_bounds, columns=columns)
    
    #prj_cols = ['ul_x', 'ul_y', 'lr_x', 'lr_y']
    #df[prj_cols] = prj_bounds
    
    return df


def split_train_test(df, pct_train=.8):
    '''
    Return two dataframes: one of training samples and the other testing. 
    '''
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
        #df_test = df[~df['row'].isin(rows) & ~df['col'].isin(cols)]

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
    trn_list = []
    oob_list = []
    keep_sets = []
    for i, df in obs:
        if df_sets.ix[i, 'n_samples'] < min_obs:
            continue
        inds = random.sample(df.index, int(len(df) * .63))
        df_b = df.ix[inds]
        df_b['set_id'] = i
        trn_list.append(df_b)
        keep_sets.append(i)
        # Get the out of bag samples
        df_oob = df[~df.index.isin(inds)]
        df_oob['set_id'] = i
        oob_list.append(df_oob)
        
    # Drop any sets that don't contain enough training observations
    #train = [df for df in dfs if len(df) > min_obs]
    df_train = pd.concat(trn_list)
    total_sets = len(df_sets)
    #keep_sets = df_train.set_id.unique()
    #keep_sets = np.array(keep_sets)
    df_drop = df_sets[~df_sets.index.isin(keep_sets)] # Sets not in keep_sets
    n_dropped = len(df_drop)
    print '%s of %s sets dropped because they contained too few observations\n' % (n_dropped, total_sets)
    df_sets = df_sets.ix[keep_sets]
    
    df_oob = pd.concat(oob_list)
    df_oob = df_oob[df_oob.set_id.isin(keep_sets)]
    
    return df_train, df_sets, df_oob, df_drop, df_test
    

def coords_to_shp(df, prj_shp, out_shp):
    '''
    Write a shapefile of rectangles from coordinates in df. Each row in df
    represents a unique set of coordinates of a rectangle feature.
    '''
    # Get spatial reference
    #import pdb; pdb.set_trace()
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

 
def get_coords(shp):
    ''' 
    Return a list of lists of the projected coordinates of vertices of shp. 
    Each list represents a feature in the dataset.
    '''    
    ds = ogr.Open(shp)
    if ds == None:
        print 'Shapefile does not exist or is not valid:\n%s' % shp
        return None
    lyr = ds.GetLayer()
    
    # For each feature, get geometry ref.
    coords = []
    for i in xrange(lyr.GetFeatureCount()): 
        this_lyr = lyr.GetFeature(i)
        geom = this_lyr.GetGeometryRef()
    
        # For each geometry ref, get coordinates from vertices
        for j in xrange(geom.GetGeometryCount()):
            this_g = geom.GetGeometryRef(j) #If geom is a Multipolgyon
            wkt = this_g.ExportToWkt()
            pts_list = wkt.replace('POLYGON','').replace('LINEARRING','').replace('(','').replace(')','').strip().split(',')
            x = [float(p.split(' ')[0]) for p in pts_list]
            y = [float(p.split(' ')[1]) for p in pts_list]
            pts = zip(x,y)
            coords.append(pts)
    
    # Need the bounds for calculating relative xy in image coords
    bounds = lyr.GetExtent()
    ds = None

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
    # Plot the support sets
    #for ind, r in df_tiles.iterrows():
        #print ind
    '''ll_x = df_tiles['ul_x']# - x_min) * x_scale
    ll_y = df_tiles['lr_y']# - y_min) * y_scale
    w = tile_size[1]# * x_scale
    h = tile_size[0]# * y_scale
    #if ll_y == minmin_y:
    sub.add_patch(plt.Rectangle((ll_x, ll_y), w, h, facecolor='none', ec='k', lw='0.5', alpha=.7, label=ind))'''
    if out_dir:
        plt.savefig(os.path.join(out_dir, 'support_sets.png'))
    else:
        plt.show()


def get_gsrd(extent_ras, cell_size, support_size, n_sets, df_train, min_obs, target_col, predict_cols, out_txt=None, shp=None, pct_train=None):
    
    print 'Generating GSRD grid...%s\n' % time.ctime(time.time())
    min_x, min_y, max_x, max_y, x_res, y_res, tx = get_raster_bounds(extent_ras)
    #ul_xy = tx[0], tx[3]
    cells = generate_gsrd_grid(cell_size, min_x, min_y, max_x, max_y, x_res, y_res)
    df_grid = pd.DataFrame(cells, columns=['ul_x', 'ul_y', 'lr_x', 'lr_y'])
    out_dir = os.path.dirname(out_txt)
    grid_shp = os.path.join(out_dir, 'gsrd_grid.shp')                         
    coords_to_shp(df_grid, shp, grid_shp)
    print 'Shapefile of support of grid cells written to:\n%s\n' % grid_shp
    
    # Get support sets
    y_size = int(support_size[0]/x_res) * x_res
    x_size = int(support_size[1]/y_res) * y_res
    support_sets = [sample_gsrd_cell(n_sets, bounds, x_size, y_size,
                                     x_res, y_res, tx) for bounds in cells]
    df_sets = pd.concat(support_sets)
    df_sets.index = xrange(len(df_sets))# Original index is range(0,n_sets)
    
    print 'Sampling observations within support sets...'
    t1 = time.time()
    #if pct_train: 
        #df_train, df_test = split_train_test(df_train, pct_train)
    df_train, df_sets, df_oob, df_drop, df_test = get_obs_within_sets(df_train, df_sets, min_obs, pct_train)
        
    set_shp = os.path.join(out_dir, 'gsrd_sets.shp')
    coords_to_shp(df_sets, shp, set_shp)
    coords_to_shp(df_drop, shp, set_shp.replace('_sets.shp', '_dropped.shp'))
    print 'Shapefile of support sets written to:\n%s' % set_shp
    print 'Time for sampling: %.1f minutes\n' % ((time.time() - t1)/60)
    
    # Plot the support sets
    '''if shp:
        print 'Plotting support sets...%s\n' % time.ctime(time.time())
        if out_txt: out_dir = os.path.join(os.path.dirname(out_txt))
        else: out_dir = None         
        coords, extent = get_coords(shp)
        set_inds = df_train.set_id.unique()
        plot_sets_on_shp(coords, 900, df_sets.ix[set_inds], support_size, out_dir)'''
    
    # Write train and test dfs to text files
    print 'Writing dataframes to disk...'
    t1 = time.time()
    train_txt = os.path.join(out_txt.replace('.txt', '_train.txt'))
    #test_txt = train_txt.replace('train', 'test')
    df_train.to_csv(train_txt, sep='\t', index=False)
    df_oob.to_csv(out_txt.replace('.txt', '_oob.txt'), sep='\t')
    if len(df_test) > 0:
        df_test.to_csv(out_txt.replace('.txt', '_test.txt'), sep='\t')
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    print 'Train and test dfs written to:\n', os.path.dirname(out_txt), '\n'
    return df_train, df_sets, df_oob


def find_file(basepath, search_str, tsa_str=None, path_filter=None):
    '''
    Return the full path within the directory tree /basepath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    '''if not os.path.exists(basepath):
        print 'basepath does not exist: \n%s' % basepath
        return None'''
     
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
    
    #import pdb; pdb.set_trace()
    if len(paths) < 1:
        #pdb.set_trace()
        sys.exit(('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str))
    
    return paths[0]


def fit_tree_classifier(x_train, y_train, max_features=None):
    ''' '''
    if not max_features: max_features=None
    dt = tree.DecisionTreeClassifier(max_features=max_features)
    dt.fit(x_train, y_train)
    
    return dt
    
    
def fit_tree_regressor(x_train, y_train, max_features=None):
    ''' '''
    if not max_features: max_features=None
    dt = tree.DecisionTreeRegressor(max_features=max_features)
    dt.fit(x_train, y_train)
    
    return dt


def fit_bdt_tree_regressor(x_samples, y_samples, pct_bagged=.63):
    
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
        pickle.dump(dt, f)
        
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
    
    # If the model is a regressor, calculate the pierson's R squared
    if model_type.lower() == 'regressor':
        oob_rate = metrics.r2_score(oob_samples, oob_prediction)
    
    # Otherwise, get the percent correct
    else:
        n_correct = len(oob_samples[oob_samples == oob_prediction])
        n_samples = len(oob_samples)
        oob_rate = int(round(float(n_correct)/n_samples * 100, 0))
    
    return oob_rate


def get_oob_rates(df_sets, df_oob, target_col, predict_cols, min_oob=0):
    
    for i, row in df_sets.iterrows():
        #with open(row.dt_file) as f:
        #    dt = pickle.load(f)
        dt = row.dt_model
        this_oob = df_oob[df_oob.set_id == i]
        oob_samples = this_oob[target_col]
        oob_predictors = this_oob[predict_cols]
        oob_rate = calc_oob_rate(dt, oob_samples, oob_predictors)
        df_sets.ix[i, 'oob_rate'] = oob_rate
    low_oob = df_sets[df_sets.oob_rate < min_oob]
    
    return df_sets, low_oob


def oob_map(ysize, xsize, nodata, nodata_mask, n_tiles, tx, support_size, df_oob, df_sets, val_col, var_cols, out_dir, file_stamp, prj, driver):
    
    t0 = time.time()
    ar_oob = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_cnt = np.full((ysize, xsize), nodata, dtype=np.int16)
    
    df_tiles, df_tiles_rc, tile_size = get_tiles(n_tiles, xsize, ysize, tx)
    total_tiles = len(df_tiles)
    
    # Find the tiles that have only nodata values
    t1 = time.time()
    print '\nFinding empty tiles...'
    empty_tiles = find_empty_tiles(df_tiles, nodata_mask, tx)
    print '%s empty tiles found of %s total tiles\n' %\
    (len(empty_tiles), total_tiles)
    # Select only tiles that are not empty
    df_tiles = df_tiles.select(lambda x: x not in empty_tiles)
    total_tiles = len(df_tiles)
    
    if 'oob_rate' not in df_sets.columns:
        for i, row in df_sets.iterrows():
            with open(row.dt_file) as f:
                dt = pickle.load(f)
            this_oob = df_oob[df_oob.set_id == i]
            oob_samples = this_oob[val_col]
            oob_predictors = this_oob[var_cols]
            df_sets.ix[i, 'oob_rate'] = calc_oob_rate(dt, oob_samples, oob_predictors)#'''
    
    for i, (t_ind, t_row) in enumerate(df_tiles.iterrows()):
        t1 = time.time()
        print 'Aggregating for %s of %s tiles' % (i + 1, total_tiles)
        print 'Tile index: ', t_ind
        
        # Calculate the size of this tile in case it's at the edge where the
        #   tile size will be slightly different
        this_size = abs(t_row.lr_y - t_row.ul_y), abs(t_row.lr_x - t_row.ul_x)
        df_these_sets = get_overlapping_sets(df_sets, t_row, this_size, support_size)
        
        rc = df_tiles_rc.ix[t_ind]
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[t_ind]
        tile_mask = nodata_mask[ul_r : lr_r, ul_c : lr_c]
        this_size = rc.lr_r - rc.ul_r, rc.lr_c - rc.ul_c
        n_sets = len(df_these_sets)
        tile_ul = t_row[['ul_x','ul_y']]
        
        # If there are no overalpping sets for this tile, fill the tile with nodata
        if t_ind in empty_tiles:
            print 'No data values for this tile...'
            continue
            this_oob = np.full(this_size, nodata, dtype=np.int16)
            this_cnt = np.full(this_size, nodata, dtype=np.int16)
        
        # Otherwise, aggregate all overlapping sets for each pixel
        else:
            print n_sets, ' Overlapping sets'
            t2 = time.time()
            
            oob_rates = []
            oob_bands = []
            for s_ind, s_row in df_these_sets.iterrows():
                #this_oob = df_oob[df_oob.set_id == s_ind]
                set_coords = s_row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
                # Fill a band for this array
                offset = calc_offset(tile_ul, (s_row.ul_x, s_row.ul_y), tx)
                #xsize = abs((s_row.ul_x - s_row.lr_x)/tx[1])
                #ysize = abs((s_row.ul_y - s_row.lr_y)/tx[5])
                tile_inds, a_inds = mosaic.get_offset_array_indices(tile_size, (support_size[0]/30, support_size[1]/30), offset)
                nrows = a_inds[1] - a_inds[0]
                ncols = a_inds[3] - a_inds[2]
                #if nrows < 0 or ncols < 0:
                #import pdb; pdb.set_trace()
                ar_oob_rate = np.full((nrows, ncols), s_row.oob_rate, dtype=np.int16)
                
                #oob_band = aggr.fill_tile_band(this_size, ar_oob_rate, tile_inds, nodata)
                oob_band = np.full(this_size, nodata)
                oob_band[tile_inds[0]:tile_inds[1], tile_inds[2]:tile_inds[3]] = ar_oob_rate
                oob_band = oob_band.astype(np.float16)
                oob_band[~tile_mask | (oob_band==nodata)] = np.nan
                oob_bands.append(oob_band)
                oob_rates.append(s_row.oob_rate)
                '''except Exception as e:
                    raise e
                    continue'''
                
            print 'Average OOB: ', int(np.mean(oob_rates))
            print 'Filling tiles: %.1f seconds' % ((time.time() - t2))
            ar_tile = np.dstack(oob_bands)
            del oob_bands
            t3 = time.time()
            this_oob = np.nanmean(ar_tile, axis=2).astype(np.int16)
            this_cnt = np.sum(~np.isnan(ar_tile), axis=2).astype(np.int16)
            
            nans = np.isnan(this_oob)
            this_oob[nans] = nodata
            this_cnt[nans] = nodata
            
            print 'Aggregating: %.1f minutes' % ((time.time() - t3)/60)
            
        print 'Total aggregation time: %.1f minutes\n' % ((time.time() - t1)/60)    
        
        # Fill in the tile in the final arrays
        #ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[t_ind]
        ar_oob[ul_r : lr_r, ul_c : lr_c] = this_oob
        ar_cnt[ul_r : lr_r, ul_c : lr_c] = this_cnt
    
    # Write final rasters to disk    
    out_path = os.path.join(out_dir, file_stamp + '_oob.bsq')
    mosaic.array_to_raster(ar_oob, tx, prj, driver, out_path, GDT_Int16, 0) 
    
    out_path = os.path.join(out_dir, out_path.replace('_oob.bsq', '_count.bsq'))
    mosaic.array_to_raster(ar_cnt, tx, prj, driver, out_path, GDT_Int16, 0)
    
    
    print '\nTotal aggregation run time: %.1f hours' % ((time.time() - t0)/3600)
    
    return ar_oob, ar_cnt, df_sets


def get_predict_array(args):

    ar = mosaic.get_mosaic(*args[1:])
    return args[0], ar.ravel()
    
    
def get_predictors(df_var, mosaic_tx, tsa_strs, tsa_ar, ar_coords, nodata_mask, out_nodata, constant_vars=None):
    '''
    Return an array of flattened predictor arrays where each predictor is a 
    separate column
    '''
    t0 = time.time()
    predictors = []
    for ind, var in enumerate(df_var.index):
        this_tsa_ar = np.copy(tsa_ar)
        data_band, search_str, basepath, by_tsa, path_filter = df_var.ix[var]
        if constant_vars: search_str = search_str.format(constant_vars['YEAR'])
        
        if by_tsa:
            files = [find_file(basepath, search_str, tsa, path_filter) for tsa in tsa_strs]
            ar_var = mosaic.get_mosaic(mosaic_tx, tsa_strs, this_tsa_ar, ar_coords, data_band, files)
            
        else:
            #this_file = find_file(basepath, search_str, path_filter=path_filter)
            #tx_, ar_var, roff_, coff_ = mosaic.get_array(this_file, data_band, ar_coords)
            try:
                this_file = find_file(basepath, search_str, path_filter=path_filter)
                tx_, ar_var, roff_, coff_ = mosaic.get_array(this_file, data_band, ar_coords)
            except:
                import pdb; pdb.set_trace()
        ar_var[nodata_mask] = out_nodata
        predictors.append(ar_var.ravel())

    #df_predict = pd.DataFrame(predictors)
    if constant_vars:
        size = predictors[0].size
        for const in sorted(constant_vars.keys()):
            val = constant_vars[const]
            predictors.append(np.full(size, val, dtype=np.int16))
    ar = np.vstack(predictors).T
    del predictors
    
    print 'Time to get predictor arrays: %.1f minutes' % ((time.time() - t0 )/60)       
    return ar
    
    
"""def predict_set(set_id, df_var, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, dt, nodata, save_stuff=None):
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
    ar_predict = get_predictors(df_var, mosaic_tx, tsa_strs, tsa_ar, ar_coords)
    del tsa_ar #Release resources from the tsa array
    
    t0 = time.time()
    nodata_mask = np.any(~(ar_predict==nodata), axis=1)
    ''' I've tried to predict in parallel here but it doesn't speed things up'''
    #p = Pool(40)
    #in_pieces = np.array_split(ar_predict[nodata_mask], 40)
    #out_pieces = p.map(par_predict, [(dt, chunk) for chunk in in_pieces])
    #import pdb; pdb.set_trace()
    #predictions = np.concatenate(out_pieces)
    predictions = dt.predict(ar_predict[nodata_mask]).astype(np.int16)
    ar_prediction = np.full(ar_predict.shape[0], nodata, dtype=np.int16)
    ar_prediction[nodata_mask] = predictions
    
    print 'Finished predicting...', time.time() - t0
    
    return ar_prediction.reshape(array_shape)"""

def predict_set(set_id, df_var, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, dt, nodata, dtype=np.int16, constant_vars=None):
    '''
    Return a predicted array for set with id==set_id
    '''
    # Get an array of tsa_ids within the bounds of ar_coords
    tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, 1, ar_coords, mosaic_tx,
                                            xsize, ysize, nodata=nodata)
    tsa_mask = tsa_ar == 0
    tsa_ar[tsa_mask] = nodata
    
    # Get the ids of TSAs this kernel covers
    tsa_ids = np.unique(tsa_ar)
    tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=nodata]
    array_shape = tsa_ar.shape

    # Get an array of predictors where each column is a flattened 2D array of a
    #   single predictor variable
    temp_nodata = -9999
    ar_predict = get_predictors(df_var, mosaic_tx, tsa_strs, tsa_ar, ar_coords, tsa_mask, temp_nodata, constant_vars)
    del tsa_ar #Release resources from the tsa array
    
    t0 = time.time()
    nodata_mask = ~ np.any(ar_predict==temp_nodata, axis=1)
    ''' I've tried to predict in parallel here but it doesn't speed things up'''
    '''p = Pool(20)
    in_pieces = np.array_split(ar_predict[nodata_mask], 20)
    out_pieces = p.map(par_predict, [(dt, chunk) for chunk in in_pieces])
    #import pdb; pdb.set_trace()
    predictions = np.concatenate(out_pieces)'''

    predictions = dt.predict(ar_predict[nodata_mask]).astype(dtype)
    ar_prediction = np.full(ar_predict.shape[0], nodata, dtype=dtype)
    ar_prediction[nodata_mask] = predictions
    
    print 'Time for predicting: %.1f minutes' % ((time.time() - t0)/60)
    
    return ar_prediction.reshape(array_shape)
    
    
def predict_set_in_pieces(set_id, df_var, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, dt, nodata, dtype=np.int16, n_pieces=10):
    '''
    Return a predicted array for set with id, set_ind
    '''
    # Get an array of tsa_ids within the bounds of ar_coords
    tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, 1, ar_coords, mosaic_tx,
                                            xsize, ysize, nodata=nodata)
    array_shape = tsa_ar.shape
    tsa_pieces = np.array_split(tsa_ar, n_pieces)
    del tsa_ar 
    
    y_res = mosaic_tx[5]
    ul_x, ul_y, lr_x, lr_y = ar_coords
    piece_coords = []
    this_uly = ul_y
    predict_pieces = []
    for i, tsa_ar in enumerate(tsa_pieces):
        print 'Predicting for piece %s of %s...' % (i + 1, n_pieces)
        t1 = time.time()
        tsa_mask = tsa_ar == 0
        tsa_ar[tsa_mask] = nodata
        
        #recalc ar_coords
        '''this_ysize, this_xsize = tsa_ar.shape
        this_uly = uly_coords[i]
        this_lry = this_uly + (this_ysize * y_res)
        these_coords = ul_x, this_uly, lr_x, this_lry'''
        
        
        #import pdb; pdb.set_trace()
        # Get the ids of TSAs this kernel covers
        tsa_ids = np.unique(tsa_ar)
        tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=nodata]
    
        # The number of rows per piece will likely be unequal, so figure out the
        #   upper and lower y for each one
        this_ysize = tsa_ar.shape[0]
        this_lry = this_uly + (this_ysize * y_res)
        these_coords = ul_x, this_uly, lr_x, this_lry
        this_uly += this_ysize * y_res        
        
        # Get an array of predictors where each column is a flattened 2D array of a
        #   single predictor variable
        ar_predict = get_predictors(df_var, mosaic_tx, tsa_strs, tsa_ar, these_coords, tsa_mask, nodata)
        
        t0 = time.time()
        nodata_mask = np.any(~(ar_predict==nodata), axis=1)
        ''' I've tried to predict in parallel here but it doesn't speed things up'''
        #p = Pool(40)
        #in_pieces = np.array_split(ar_predict[nodata_mask], 40)
        #out_pieces = p.map(par_predict, [(dt, chunk) for chunk in in_pieces])
        #import pdb; pdb.set_trace()
        #predictions = np.concatenate(out_pieces)
    
        predictions = dt.predict(ar_predict[nodata_mask]).astype(dtype)
        ar_prediction = np.full(ar_predict.shape[0], nodata, dtype=dtype)
        ar_prediction[nodata_mask] = predictions
        predict_pieces.append(ar_prediction)
        print 'Time for predicting this piece: %.1f minutes' % ((time.time() - t1)/60)
    try:
        ar_prediction = np.concatenate(predict_pieces)
    except:
        import pdb; pdb.set_trace()
    print 'Time for predicting: %.1f minutes' % ((time.time() - t0)/60)
    
    return ar_prediction.reshape(array_shape)


def par_predict(args):
    '''Helper function to parallelize predicting for a decision tree'''
    dt, ar = args
    try:
        return dt.predict(ar)
    except:
        sys.exit(traceback.print_exception(*sys.exc_info()))
        

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
        ar_pred = ar_pred.astype(np.float32)
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
    t3 = time.time()
    ar_out = np.full(ar.shape, nodata, dtype=np.int16)
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
    
    return np.round(w_mean,0).astype(np.int16)
    

def pct_vote(ar, ar_vote, ar_count):
    
    ar_eq = ar == ar_vote
    ar_sum = ar_eq.sum(axis=2)
    ar_pct = np.round(ar_sum/ar_count.astype(np.float16) * 100).astype(np.uint8)
    
    return ar_pct
    

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
    Set ar == mask_val where mask is true
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
    ar_mean = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_vote = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_stdv = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_coun = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_impr = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_wtmn_10 = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_wtmn_20 = np.full((ysize, xsize), nodata, dtype=np.int16)#'''
    
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
            this_mean = np.full(this_size, nodata, dtype=np.int16)
            this_vote = np.full(this_size, nodata, dtype=np.int16)
            this_stdv = np.full(this_size, nodata, dtype=np.int16)
            this_coun = np.full(this_size, nodata, dtype=np.int16)
            this_impr = np.full(this_size, nodata, dtype=np.int16)
            this_wtmn_10 = np.full(this_size, nodata, dtype=np.int16)
            this_wtmn_20 = np.full(this_size, nodata, dtype=np.int16)
        
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
                ar_import = np.full(ar_pred.shape, max_import, dtype=np.int16)
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
        ar_mean[ul_r : lr_r, ul_c : lr_c] = this_mean.astype(np.int16)
        ar_vote[ul_r : lr_r, ul_c : lr_c] = this_vote.astype(np.int16)
        ar_stdv[ul_r : lr_r, ul_c : lr_c] = this_stdv.astype(np.int16)
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
    
"""def aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_ds, support_size, prediction_dir, df_sets, out_dir, file_stamp, prj, driver, mosaic_nodata=0):
    
    t0 = time.time()
    ar_mean = np.full((ysize, xsize), nodata, dtype=np.int16)
    #ar_vote = np.full((ysize, xsize), nodata, dtype=np.int16)
    #ar_stdv = np.full((ysize, xsize), nodata, dtype=np.int16)
    #ar_coun = np.full((ysize, xsize), nodata, dtype=np.int16)
    #ar_impr = np.full((ysize, xsize), nodata, dtype=np.int16)
    #ar_wtmn_10 = np.full((ysize, xsize), nodata, dtype=np.int16)#'''
    #ar_wtmn_20 = np.full((ysize, xsize), nodata, dtype=np.int16)
    
    #mosaic_ds = gdal.Open(mosaic_path)
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
    '''t1 = time.time()
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
            this_mean = np.full(this_size, nodata, dtype=np.int16)
            #this_vote = np.full(this_size, nodata, dtype=np.int16)
            #this_stdv = np.full(this_size, nodata, dtype=np.int16)
            #this_coun = np.full(this_size, nodata, dtype=np.int16)
            #this_impr = np.full(this_size, nodata, dtype=np.int16)
            #this_wtmn_10 = np.full(this_size, nodata, dtype=np.int16)
            #this_wtmn_20 = np.full(this_size, nodata, dtype=np.int16)
        
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
                '''try:
                    with open(s_row.dt_file, 'rb') as f: 
                        dt_model = pickle.load(f)
                    #max_importance, importance = get_max_importance(dt_model)
                    #importance_list.append(importance)
                    ar_import = np.full(ar_pred.shape, s_row.max_importance, dtype=np.int16)
                    #max_importance = important_features(dt_model, ar_pred, nodata)
                    import_band = fill_tile_band(this_size, ar_import, tile_inds, nodata)
                    importance_bands.append(import_band)
                except Exception as e:
                    print e
                    continue#'''
            print 'Filling tiles: %.1f seconds' % ((time.time() - t2))
            #import pdb; pdb.set_trace()
            ar_tile = np.dstack(pred_bands)
            #nd_impr = np.dstack(importance_bands)
            del pred_bands#, importance_bands#, ar_import
            t3 = time.time()
            this_mean = np.nanmean(ar_tile, axis=2)
            #this_vote = mode(ar_tile, axis=2)
            #this_stdv = np.nanstd(ar_tile, axis=2) * 100 #Multiply b/c converting to int
            #this_coun = np.sum(~np.isnan(ar_tile), axis=2)
            #this_impr = mode(nd_impr, axis=2)
            #this_wtmn_10 = weighted_mean(ar_tile, this_vote, c=10)
            #this_wtmn_20 = weighted_mean(ar_tile, this_vote, c=20)
            
            nans = np.isnan(this_mean)
            nans = np.isnan(this_vote)
            this_mean[nans] = nodata
            #this_stdv[nans] = nodata
            #this_impr[nans] = nodata
            #this_vote[nans] = nodata
            #this_coun[nans] = nodata
            #this_wtmn_10[nans] = nodata
            #this_wtmn_20[nans] = nodata
            
            print 'Aggregating: %.1f minutes' % ((time.time() - t3)/60)
        print 'Total aggregation time: %.1f minutes\n' % ((time.time() - t1)/60)    
        
        # Fill in the tile in the final arrays
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[t_ind]
        ar_mean[ul_r : lr_r, ul_c : lr_c] = this_mean.astype(np.int16)
        #ar_vote[ul_r : lr_r, ul_c : lr_c] = this_vote.astype(np.int16)
        #ar_stdv[ul_r : lr_r, ul_c : lr_c] = this_stdv.astype(np.int16)
        #ar_coun[ul_r : lr_r, ul_c : lr_c] = this_coun
        #ar_impr[ul_r : lr_r, ul_c : lr_c] = this_impr
        #ar_wtmn_10[ul_r : lr_r, ul_c : lr_c] = this_wtmn_10
        #ar_wtmn_20[ul_r : lr_r, ul_c : lr_c] = this_wtmn_20
    
    # Mask arrays
    #mask, tx_mask = get_nonforest_mask(lc_path, ag_path, lc_vals, ag_vals)
    #mask_array(ar_mean, mask, mosaic_tx, tx_mask)
    #mask_array(ar_vote, mask, mosaic_tx, tx_mask)
    
    # Write final rasters to disk
    out_template = os.path.join(out_dir, file_stamp + '_%s.bsq')
    out_path = out_template % 'mean'
    mosaic.array_to_raster(ar_mean, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata)
    
    out_path = out_template % 'vote'
    #mosaic.array_to_raster(ar_vote, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata)   
    
    out_path = out_template % 'stdv'
    #mosaic.array_to_raster(ar_stdv, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata)
    
    #out_path = out_path.replace('stdv', 'countagg')
    #mosaic.array_to_raster(ar_coun, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = out_template % 'importance'
    #mosaic.array_to_raster(ar_impr, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata)  
    
    #out_path = os.path.join(out_dir, file_stamp + '_weightedmean_10.bsq')
    #mosaic.array_to_raster(ar_wtmn_10, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = os.path.join(out_dir, file_stamp + '_weightedmean_20.bsq')
    #mosaic.array_to_raster(ar_wtmn_20, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    print '\nTotal aggregation run time: %.1f hours' % ((time.time() - t0)/3600)
    #return predictions, df_sets, df_train
    #return df_tiles, df_these_sets, ar_mean, ar_tile, ar_out
    #return ar_out
    #del ar_impr, ar_pred, ar_stdv, ar_tile, this_impr, this_mean, this_stdv, this_vote, #this_wtmn_20, ar_coun, this_coun, this_wtmn_10, ar_wtmn_20, ar_wtmn_10
    
    return None, ar_vote, None, None#"""


def aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_ds, support_size, prediction_dir, df_sets, out_dir, file_stamp, prj, driver, mosaic_nodata=0):
    
    t0 = time.time()
    #ar_mean = np.full((ysize, xsize), nodata, dtype=np.int16)
    ar_vote = np.full((ysize, xsize), nodata, dtype=np.uint8)
    #ar_stdv = np.full((ysize, xsize), nodata, dtype=np.int16)
    #ar_coun = np.full((ysize, xsize), nodata, dtype=np.uint8)
    ar_impr = np.full((ysize, xsize), nodata, dtype=np.uint8)
    ar_pcvt = np.full((ysize, xsize), nodata, dtype=np.uint8)
    #ar_wtmn_10 = np.full((ysize, xsize), nodata, dtype=np.int16)#'''
    #ar_wtmn_20 = np.full((ysize, xsize), nodata, dtype=np.int16)
    
    #mosaic_ds = gdal.Open(mosaic_path)
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
        
        rc = df_tiles_rc.ix[t_ind]
        this_size = rc.lr_r - rc.ul_r, rc.lr_c - rc.ul_c
        n_sets = len(df_these_sets)
        
        # Load overlapping predictions from disk and read them as arrays
        tile_ul = t_row[['ul_x','ul_y']]
        predictions = load_predictions(prediction_dir, df_these_sets, tile_ul, this_size)
        
        # If there are no overalpping sets for this tile, fill the tile with nodata
        '''if t_ind in empty_tiles:
            print 'No overlapping sets for this tile'
            continue
            #this_mean = np.full(this_size, nodata, dtype=np.int16)
            this_vote = np.full(this_size, nodata, dtype=np.int16)
            #this_stdv = np.full(this_size, nodata, dtype=np.int16)
            this_coun = np.full(this_size, nodata, dtype=np.int16)
            this_impr = np.full(this_size, nodata, dtype=np.int16)
            #this_wtmn_10 = np.full(this_size, nodata, dtype=np.int16)
            #this_wtmn_20 = np.full(this_size, nodata, dtype=np.int16)'''
        
        # Otherwise, aggregate all overlapping sets for each pixel
        #else:
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
                ar_import = np.full(ar_pred.shape, s_row.max_importance, dtype=np.uint8)
                #max_importance = important_features(dt_model, ar_pred, nodata)
                import_band = fill_tile_band(this_size, ar_import, tile_inds, nodata)
                importance_bands.append(import_band)
            except Exception as e:
                print e
                continue#'''
            print 'Filling tiles: %.1f seconds' % ((time.time() - t2))
            
            ar_tile = np.dstack(pred_bands)
            nd_impr = np.dstack(importance_bands)
            del pred_bands, importance_bands#, ar_import
            t3 = time.time()
            #this_mean = np.nanmean(ar_tile, axis=2)
            this_vote = mode(ar_tile, axis=2)
            #this_stdv = np.nanstd(ar_tile, axis=2) * 100 #Multiply b/c converting to int
            this_coun = np.sum(~np.isnan(ar_tile), axis=2)
            this_impr = mode(nd_impr, axis=2)
            this_pcvt = pct_vote(ar_tile, this_vote, this_coun)
            #this_wtmn_10 = weighted_mean(ar_tile, this_vote, c=10)
            #this_wtmn_20 = weighted_mean(ar_tile, this_vote, c=20)
            
            nans = np.isnan(this_vote)
            #this_mean[nans] = nodata
            #this_stdv[nans] = nodata
            this_impr[nans] = nodata
            this_vote[nans] = nodata
            #this_coun[nans] = nodata
            this_pcvt[nans] = nodata
            #this_wtmn_10[nans] = nodata
            #this_wtmn_20[nans] = nodata
            
            print 'Aggregating: %.1f minutes' % ((time.time() - t3)/60)
        print 'Total aggregation time: %.1f minutes\n' % ((time.time() - t1)/60)    
        
        # Fill in the tile in the final arrays
        ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[t_ind]
        #ar_mean[ul_r : lr_r, ul_c : lr_c] = this_mean.astype(np.int16)
        ar_vote[ul_r : lr_r, ul_c : lr_c] = this_vote.astype(np.uint8)
        #ar_stdv[ul_r : lr_r, ul_c : lr_c] = this_stdv.astype(np.int16)
        #ar_coun[ul_r : lr_r, ul_c : lr_c] = this_coun
        ar_impr[ul_r : lr_r, ul_c : lr_c] = this_impr
        ar_pcvt[ul_r : lr_r, ul_c : lr_c] = this_pcvt
        #ar_wtmn_10[ul_r : lr_r, ul_c : lr_c] = this_wtmn_10
        #ar_wtmn_20[ul_r : lr_r, ul_c : lr_c] = this_wtmn_20
    
    #Mask arrays
    #ar_mean[~mask] = nodata
    ar_vote[~mask] = nodata
    ar_impr[~mask] = nodata
    #ar_coun[~mask] = nodata
    ar_pcvt[~mask] = nodata
    #ar_stdv[~mask] = nodata * 100
    
    # Write final rasters to disk
    out_template = os.path.join(out_dir, file_stamp + '_%s.bsq')
    out_path = out_template % 'mean'
    #mosaic.array_to_raster(ar_mean, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata)
    
    out_path = out_template % 'vote'
    mosaic.array_to_raster(ar_vote, mosaic_tx, prj, driver, out_path, GDT_Byte, nodata)   
    
    out_path = out_template % 'stdv'
    #mosaic.array_to_raster(ar_stdv, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata * 100)
    
    #out_path = out_path.replace('stdv', 'countagg')
    #mosaic.array_to_raster(ar_coun, mosaic_tx, prj, driver, out_path, GDT_Int16, nodata)
    
    out_path = out_template % 'importance'
    mosaic.array_to_raster(ar_impr, mosaic_tx, prj, driver, out_path, GDT_Byte, nodata)  
    
    out_path = out_template % 'pct_vote'
    mosaic.array_to_raster(ar_pcvt, mosaic_tx, prj, driver, out_path, GDT_Byte, nodata) 
    
    #out_path = os.path.join(out_dir, file_stamp + '_weightedmean_10.bsq')
    #mosaic.array_to_raster(ar_wtmn_10, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    out_path = os.path.join(out_dir, file_stamp + '_weightedmean_20.bsq')
    #mosaic.array_to_raster(ar_wtmn_20, mosaic_tx, prj, driver, out_path, GDT_Int32, nodata)
    
    print '\nTotal aggregation run time: %.1f hours' % ((time.time() - t0)/3600)
    #return predictions, df_sets, df_train
    #return df_tiles, df_these_sets, ar_mean, ar_tile, ar_out
    #return ar_out
    #del ar_impr, ar_pred, ar_stdv, ar_tile, this_impr, this_mean, this_stdv, #this_vote, #this_wtmn_20, ar_coun, this_coun, this_wtmn_10, ar_wtmn_20, ar_wtmn_10
    
    ar_mean = None
    return ar_mean, ar_vote, pct_import, df_sets#"""
    

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
