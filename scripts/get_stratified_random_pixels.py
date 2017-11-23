# -*- coding: utf-8 -*-
"""
Generate stratified samples of specified bin sizes from a given input raster 

@author: Sam Hooper, samhooperstudio@gmail.com
"""

import random
import sys
import os
import time
import shutil
import pandas as pd
import numpy as np
from datetime import datetime
from osgeo import gdal, ogr

import extract_xy_by_mosaic as extract
import stem

#gdal.SetCacheMax(2 * 30)

def read_params(txt):
    '''
    Return a dictionary from parsed parameters in txt
    '''
    d = {}
    
    # Read in the text file
    try:
        with open(txt) as f:
            input_vars = [line.split(";") for line in f]
    except: 
        print 'Problem reading parameter file: ', txt
        return None
    
    # Set the dictionary key to whatever is left of the ";" and the val
    #   to whatever is to the right. Strip whitespace too.
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
                '"{0}"'.format(var[1].strip(" ").replace("\n", ""))
        
    print 'Parameters read from:\n', txt, '\n'
    return d


def parse_bins(bin_str):
    ''' Integerize a string of min:max and return as a list of length 2'''
    bin_list = [b.split(':') for b in bin_str.split(',')]
    bins = [(int(mn), int(mx)) for mn, mx in bin_list]
    
    return bins


def extract_by_kernel(ar, rows, cols, data_type, col_name, nodata):
    
    row_dirs = [-1,-1,-1, 0, 0, 0, 1, 1, 1]
    col_dirs = [-1, 0, 1,-1, 0, 1,-1, 0, 1]
    kernel_rows = [row + d + 1 for row in rows for d in row_dirs]
    kernel_cols = [col + d + 1 for col in cols for d in col_dirs]
    
    ar_buf = np.full([dim + 2 for dim in ar.shape], nodata, dtype=np.int32)
    ar_buf[1:-1, 1:-1] = ar
    del ar
    
    kernel_vals = ar_buf[kernel_rows, kernel_cols].reshape(len(rows), len(row_dirs))
    train_stats = pd.DataFrame(extract.calc_row_stats(kernel_vals, data_type, col_name, nodata))
    vals = train_stats[col_name].astype(np.int32)
    
    return vals


def calc_proportional_area(df_tiles, shp):
    
    ds = ogr.Open(shp)
    lyr = ds.GetLayer()
    feat = lyr.GetFeature(0)
    shp_boundary = feat.GetGeometryRef()
    
    for i, coords in df_tiles.iterrows():
        # Create geometry
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(coords.ul_x, coords.ul_y) #ul vertex
        ring.AddPoint(coords.lr_x, coords.ul_y) #ur vertex
        ring.AddPoint(coords.lr_x, coords.lr_y) #lr vertex
        ring.AddPoint(coords.ul_x, coords.lr_y) #ll vertex
        ring.CloseRings()
        tile = ogr.Geometry(ogr.wkbPolygon)
        tile.AddGeometry(ring)
        
        #Calulate intersection
        try:
            intersection = tile.Intersection(shp_boundary)
            area = intersection.GetArea()
            df_tiles.ix[i, 'area'] = area
        except:
            print 'Problem with tile', i
    df_tiles['pct_area'] = df_tiles.area/df_tiles.area.sum()
    
    ds = None
        
    #return df_tiles


def calc_strata_sizes(ar, nodata_mask, sample_size, bins, scheme='proportional', bin_scale=1, zero_inflation=0):
    
    n_pixels = float(ar[nodata_mask].size)

    if scheme == 'proportional':
        counts, _ = np.histogram(ar, bins=[b[1] for b in bins] + [bins[-1][1] + 1])# bin param for np.histogram is all left edge ecxept last is right edge
        bin_percents = counts.astype(np.float)/n_pixels
        pcts_center = (bin_percents.max() - bin_percents.min())/2.0 #middle of range
        scaled_dif = (bin_percents - pcts_center)/float(bin_scale)
        scaled_percents = scaled_dif + pcts_center
        scaled_percents /= scaled_percents.sum()
        # rescale all vals so sum for this tile == samples_per
        scaled_samples_per = (scaled_percents * sample_size)
        n_per_bin = scaled_samples_per * float(sample_size)/sum(scaled_samples_per)
    
    elif scheme == 'equal':
        n_bins = len(bins)
        bin_size = float(sample_size)/n_bins
        n_per_bin = np.repeat(bin_size, n_bins)
        if zero_inflation:
            bin_maxes = np.array([b[1] for b in bins])
            bin0_ind = np.arange(n_bins)[bin_maxes == 0]
            n_per_bin[bin0_ind] *= zero_inflation
        scaled_percents = n_per_bin/sum(n_per_bin)
        n_per_bin = scaled_percents * sample_size

    return n_per_bin, scaled_percents
    

def get_stratified_sample_by_tile(raster_path, col_name, data_band, n_samples, bins, min_sample=None, max_sample=None, pct_train=1, nodata=None, sampling_scheme='equal', zero_inflation=None, data_type='continuous', kernel=False, n_tiles=(1, 1), boundary_shp=None, bin_scale=1, n_per_tile=None):
    '''
    Return a dataframe of stratified randomly sampled pixels from raster_path
    '''
    print 'Reading the raster_path... %s\n' % datetime.now()
    ds = gdal.Open(raster_path)
    tx = ds.GetGeoTransform()
    band = ds.GetRasterBand(data_band)
    #ar_full = band.ReadAsArray()
    #ar_data = band.ReadAsArray()
    if nodata == None:
        nodata = band.GetNoDataValue()
        if nodata == None:
            sys.exit('Could not obtain nodata value from dataset and' +\
            ' none specified in parameters file. Try re-running with' +\
            'nodata specified.')
    
    # Buffer the array by 1 on each side to be able to extract 3x3 kernels
    #ar = np.full([dim + 2 for dim in ar_data.shape], nodata, dtype=np.int32)
    #ar[1:-1, 1:-1] = ar_data
    xsize = ds.RasterXSize 
    ysize = ds.RasterYSize
    df_tiles, df_tiles_rc, tile_size = stem.get_tiles(n_tiles, xsize, ysize, tx)
    if boundary_shp:
        boundary_ds = ogr.Open(boundary_shp)
        boundary_lyr = boundary_ds.GetLayer()
        empty_tiles = stem.find_empty_tiles(df_tiles, boundary_lyr, tx, nodata=0)
        df_tiles.drop(empty_tiles, inplace=True)
        df_tiles_rc.drop(empty_tiles, inplace=True)
        total_tiles = len(df_tiles)
        
    # If a boundary shapefile is given, calculate the proportional area within
    #   each tile
    total_tiles = len(df_tiles)
    if not n_per_tile:
        if boundary_shp:
            calc_proportional_area(df_tiles, boundary_shp) # Calcs area in place
            df_tiles_rc['n_samples'] = n_samples * df_tiles.pct_area
        else:
            df_tiles_rc['n_samples'] = n_samples/total_tiles
    else:
        df_tiles_rc['n_samples'] = n_per_tile
    
    # For each tile, get random samples for each bin
    train_rows = []
    train_cols = []
    test_rows = []
    test_cols = []
    empty_tiles = []
    for i, tile_coords in df_tiles_rc.iterrows():
        t1 = time.time()
        print 'Sampling for tile %s of %s...' % (i + 1, total_tiles)
        if tile_coords.n_samples == 0:
            print '\tSkipping this tile because all pixels == nodata...\n'
            empty_tiles.append(i)
            continue
        ul_r, lr_r, ul_c, lr_c = tile_coords[['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        tile_ysize = lr_r - ul_r
        tile_xsize = lr_c - ul_c 

        ar = band.ReadAsArray(ul_c, ul_r, tile_xsize, tile_ysize)
        if not type(ar) == np.ndarray:
            import pdb; pdb.set_trace()
        nodata_mask = ar != nodata
        if not nodata_mask.any():
            print '\tSkipping this tile because all pixels == nodata...\n'
            empty_tiles.append(i)
            continue
        ar_rows, ar_cols = np.indices(ar.shape)
        ar_rows = ar_rows + ul_r
        ar_cols = ar_cols + ul_c
        
        n_per_bin, scaled_pcts = calc_strata_sizes(ar, nodata_mask, tile_coords.n_samples, bins, sampling_scheme, bin_scale, zero_inflation)
        
        for i, (this_min, this_max) in enumerate(bins):
            #t2 = time.time()
            
            try:
                this_sample_size = n_per_bin[i]
                if min_sample:
                    this_sample_size = int(max(n_per_bin[i], min_sample))
                if max_sample:
                    this_sample_size = int(min(n_per_bin[i], max_sample))
            except:
                import pdb; pdb.set_trace()
            print 'Sampling between %s and %s: %s pixels (%.1f%% of sample for this tile) ' % (this_min, this_max, this_sample_size, scaled_pcts[i] * 100)
            mask = (ar > this_min) & (ar <= this_max) & nodata_mask
            these_rows = ar_rows[mask]
            these_cols = ar_cols[mask]
            
            '''if this_max == 0 and zero_inflation:
                this_sample_size *= zero_inflation'''
            
            # If there aren't enough pixels to generate samples for this bin
            if these_rows.size < this_sample_size:
                #import pdb; pdb.set_trace()
                print (('Not enough pixels between {0} and {1} to generate {2}' +\
                ' random samples. Returning all {3} pixels for this bin.'))\
                .format(this_min, this_max, this_sample_size, these_rows.size)
                tr_rows = these_rows
                tr_cols = these_cols
            else:
                samples = random.sample(xrange(len(these_rows)), this_sample_size)
                tr_rows = these_rows[samples]
                tr_cols = these_cols[samples]


            
            # If pct_train is specified, split the sample indices into train/test sets
            te_rows = []
            te_cols = []
            if pct_train:
                split_ind = len(tr_rows) * pct_train
                te_rows = tr_rows[split_ind:]
                te_cols = tr_cols[split_ind:]
                tr_rows = tr_rows[:split_ind]
                tr_cols = tr_cols[:split_ind]
    
            train_rows.extend(tr_rows)
            train_cols.extend(tr_cols)
            test_rows.extend(te_rows)
            test_cols.extend(te_cols)
            #print '%.1f seconds\n' % (time.time() - t2)
        print 'Time for this tile: %.1f minutes\n' % ((time.time() - t1)/60)
    del tr_rows, tr_cols, te_rows, te_cols, ar
    
    # Read the whole raster in to extract stuff
    ar = band.ReadAsArray()
    
    # If True, extract with 3x3 kernel. Otherwise, just get the vals (row,col)
    if kernel:
        train_vals = extract_by_kernel(ar, train_rows, train_cols, data_type, col_name, nodata)
    else:
        train_vals = ar[train_rows, train_cols]
    
    # Calculate x and y for later extractions
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    train_x = [int(ul_x + c * x_res) for c in train_cols]
    train_y = [int(ul_y + r * y_res) for r in train_rows]
    df_train = pd.DataFrame({'x': train_x,
                             'y': train_y,
                             'row': train_rows,
                             'col': train_cols,
                             col_name: train_vals
                             })
                             
    # If training and testing samples were split, get test vals                    
    df_test = None
    if pct_train > 1:
        if kernel:
            test_vals = extract_by_kernel(ar, train_rows, train_cols, data_type, col_name, nodata)
        else:
            test_vals = ar[test_rows, test_cols]
        test_x = [int(ul_x + c * x_res) for c in test_cols]
        test_y = [int(ul_y + r * y_res) for r in test_rows]
        df_test = pd.DataFrame({'x': test_x, 
                                'y': test_y, 
                                'row': test_rows,
                                'col': test_cols,
                                col_name: test_vals
                                })
    ds = None
    
    # In case empty tiles weren't filtered out already
    df_tiles = df_tiles.ix[~df_tiles.index.isin(empty_tiles)]
    
    return df_train, df_test, df_tiles


def main(params, data_band=1, nodata=None, sampling_scheme='proportional', data_type='continuous', kernel=False, boundary_shp=None, bin_scale=1, min_sample=None, max_sample=None, n_samples=None, n_per_tile=None):
    
    t0 = time.time()
    data_band = None
    nodata = None
    zero_inflation = None
    
    # Read params and make variables from each line
    inputs = read_params(params)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    
    #out_dir = os.path.dirname(out_txt)
    '''if not os.path.exists(out_dir):
        print 'Warning: output directory does not exist. Creating directory...'
        os.makedirs(out_dir)'''
    
    # Integerize numeric params
    if 'data_band' in locals(): data_band = int(data_band)
    if 'nodata' in locals(): nodata = int(nodata)
    if 'pct_train' in locals(): 
        pct_train = float(pct_train)
    else: 
        pct_train = None
    if zero_inflation: zero_inflation = int(zero_inflation)
    if 'bin_scale' in inputs:
        bin_scale = float(bin_scale)
    if 'min_sample' in inputs:
        min_sample = int(min_sample)
    if 'max_sample' in inputs:
        max_sample = int(max_sample)
    if 'n_per_tile' in inputs:
        n_per_tile = int(n_per_tile)
    if 'n_samples' in inputs:
        n_samples = int(n_samples)
    try:
       bins = parse_bins(bins)
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    # If number of tiles not given, need to calculate them
    if 'n_tiles' in inputs:
        n_tiles = [int(i) for i in n_tiles.split(',')]
    else:
        n_tiles = 3, 10
        print 'Using default tile size of %s x %s ....' % n_tiles
        
    
    # Generate samples
    df_train, df_test, df_tiles = get_stratified_sample_by_tile(raster_path, col_name,
                                                      data_band, n_samples,
                                                      bins, min_sample, max_sample, pct_train, nodata,
                                                      sampling_scheme, zero_inflation, data_type,
                                                      kernel, n_tiles, boundary_shp, 
                                                      bin_scale=bin_scale, n_per_tile=n_per_tile)
    df_train['obs_id'] = df_train.index
    
    # Write samples to text file
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    bn = '{0}_{1}_sample_{2}_{3}_{4}.txt'.format(col_name, sampling_scheme, len(df_train), date_str, time_str)
    #bn = os.path.basename(out_txt)
    stamp = bn[:-4]
    out_dir = os.path.join(out_dir, stamp)
    #if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    out_txt = os.path.join(out_dir, bn)
    df_train.to_csv(out_txt, sep='\t', index=False)
    print 'Sample written to:\n%s\n' % out_txt
    
    shutil.copy2(params, out_dir) #Copy the params for reference
    
    if pct_train > 1:
        df_test['obs_id'] = df_test.index
        test_txt = out_txt.replace('%s.txt' % stamp, '%s_test.txt' % stamp)
        df_test.to_csv(test_txt, sep='\t', index=False)
        print 'Test samples written to directory:\n%s' % out_dir
    
    if n_tiles != [1,1]:
        if boundary_shp:
            out_shp = os.path.join(out_dir, 'sampling_tiles.shp')
            stem.coords_to_shp(df_tiles, boundary_shp, out_shp)
        else:
            tile_txt = os.path.join(out_dir, 'sampling_tiles.txt')
            df_tiles.to_csv(tile_txt, sep='\t', index=False)
            
    print '\nTotal time: %.1f minutes' % ((time.time() - t0)/60)
    
    return out_txt
        
        
if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))

#raster_path = '/vol/v2/stem/canopy/truth_map/nlcd_canopy_2001_nodata.tif'
#out_txt = '/home/server/student/homes/shooper/scripts/random_pixel_test.txt'
#params = '/vol/v2/stem/scripts/get_stratified_random_pixels_params.txt'
#df = main(params)
