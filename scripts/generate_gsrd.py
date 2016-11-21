# -*- coding: utf-8 -*-
"""
Script to generate geographically stratified random design (GSRD) grid

-get study area bounds
-randomize position of grid
    -start with seed ul and 
-evenly and randomly sample points in each statum
-calculate bounds
-store ul_x, ul_y, lr_x, and lr_y in separate columns of a df
-randomly sample x's and y's within each set of bounds to get support set centers
-create a dataframe of support sets 

"""
import gdal
import ogr
#import osr
import os
import random
import matplotlib
import numpy as np
import pandas as pd
import time
import matplotlib
matplotlib.use('Agg')
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
#from matplotlib.patches import Polygon

# Turn off the annoying setting value on a copy warning
pd.options.mode.chained_assignment = None

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
    
    # Calculate bounding coords for from support set centers and make sure 
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
    df = pd.DataFrame(these_bounds, columns=['ul_x', 'ul_y', 'lr_x', 'lr_y', 'ctr_x', 'ctr_y'])
    
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
    
    return df_train, df_sets, df_oob, df_drop
    

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
    
    print 'Sampling observations within support sets...%s\n' % time.ctime(time.time())
    t1 = time.time()
    #if pct_train: 
        #df_train, df_test = split_train_test(df_train, pct_train)
    df_train, df_sets, df_oob, df_drop = get_obs_within_sets(df_train, df_sets, min_obs, pct_train)
        
    set_shp = os.path.join(out_dir, 'gsrd_sets.shp')
    coords_to_shp(df_sets, shp, set_shp)
    coords_to_shp(df_drop, shp, set_shp.replace('_sets.shp', '_dropped.shp'))
    print 'Shapefile of support sets written to:\n%s\n' % set_shp
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
    df_oob.to_csv(train_txt.replace('_train.txt', '_oob.txt'), sep='\t')
    #df_sets.to_csv(train_txt.replace('_train.txt', '_sets.txt'), sep='\t', index=False)
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    print 'Train and test dfs written to:\n', os.path.dirname(out_txt), '\n'
    return df_train, df_sets, df_oob

''' Testing '''
"""
extent_ras = '/vol/v1/general_files/datasets/spatial_data/CAORWA_TSA_lt_only.bsq'
#cell_size = (150000, 100000)
#support_size = (200000, 150000)
cell_size = (300000, 200000)
support_size = (400000, 300000)
sets_per_cell = 5
min_obs = 25
pct_train = 1
xy_txt = '/vol/v2/stem/canopy/samples/canopy_sample10000_20160308_2046/canopy_sample10908_20160308_2046_predictors.txt'
#xy_txt = '/vol/v2/stem/scripts/xy_prj_tsa.txt'
shp = '/vol/v2/stem/extent_shp/orwaca.shp'
out_txt = '/vol/v2/stem/imperv/models/delete/gsrd_test_delete.txt'
target_name = ''


train, sets = get_gsrd(extent_ras, cell_size, support_size, sets_per_cell, xy_txt, min_obs, target_name, out_txt, shp=shp)

#extent_ras, cell_size, support_size, sets_per_cell, xy_txt, min_obs, target_name, out_txt, shp=shp#"""
    

    
    