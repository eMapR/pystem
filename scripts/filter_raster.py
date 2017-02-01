# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:13:57 2016

@author: shooper
"""
import os
import sys
import gdal
import time
import shutil
import numpy as np
from multiprocessing import Pool
from scipy import ndimage as ndi
from stem import get_tiles, find_empty_tiles, mode, coords_to_shp
from mosaic_by_tsa import array_to_raster


def read_params(txt):
    '''
    Return a dictionary and a dataframe from parsed parameters in txt
    '''
    if not os.path.exists(txt):
        print 'Param file does not exist:\n%s' % txt
        return None

    # Read in the rest of the text file line by line
    d = {}
    try:
        with open(txt) as f:
            input_vars = [line.split(";") for line in f]
    except:
        print 'Problem reading parameter file:\n', txt
        return None

    # Set the dictionary key to whatever is left of the ";" and the val
    #   to whatever is to the right. Strip whitespace too.
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] = var[1].strip(' ').replace('\n', '')

    print '\nParameters read from:\n', txt, '\n'
    
    return d


def circle_mask(diameter):
    
    radius = diameter/2
    kernel = np.zeros((2 * radius + 1, 2 * radius + 1))
    y, x = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    mask = x**2 + y**2 <= radius**2    
    kernel[mask] = 1

    return kernel


def pct_nonzero(ar):
    
    mask = ~np.isnan(ar)
    ar = ar[mask]
    n_pixels = float(ar.size)
    if n_pixels == 0:
        return 0
    pct = np.count_nonzero(ar)/n_pixels * 100
    
    return pct
    

def is_equal_to(ar, center_idx):
    
    center_val = ar[center_idx]
    ar = ar[~np.isnan(ar)]
    if np.all(ar == center_val):
        return 1
    else:
        return 0


def par_filter(args):
    
    t0 = time.time()
    ind, ar, func, kernel, extra_args, i, n_tiles = args
    ar = ndi.generic_filter(ar, func, footprint=kernel, extra_arguments=extra_args)
    print 'Time for tile %s of %s: %.1f minutes' % (i, n_tiles, ((time.time() - t0)/60))

    return ind, ar


def main(params, n_tiles=(25, 15), n_jobs=20, kernel_type='circle', filter_value=None):
    
    t0 = time.time()
    
    # Read params and make variables from text
    inputs = read_params(params)
        
    # Check params
    try:
        path = inputs['path']
        function = inputs['function']
        out_path = inputs['out_path']
        kernel_size = int(inputs['kernel_size'])
        databand = int(inputs['databand'])
    except KeyError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    if 'n_jobs' in inputs: n_jobs = int(inputs['n_jobs'])
    if 'n_tiles' in inputs: n_tiles = [int(n) for n in inputs['n_tiles'].split(',')]
    if 'nodata' in inputs: nodata = int(inputs['nodata'])
    
    extra_args = () # The default for ndi.generic_filter 'extra_args' is an empty tuple
    if 'average' in function.lower():
        func = np.nanmean
    elif 'mode' in function.lower():
        func = mode
    elif 'area' in function.lower():
        func = pct_nonzero
        if not filter_value and not 'filter_value' in inputs:
            sys.exit('Cannot calculate percent area without filter_value. ' +\
            'Try specifying filter_value in parameters file.')
        else:
            filter_value = int(inputs['filter_value'])
    elif 'equal' in function.lower():
        func = is_equal_to
        center_idx = kernel_size**2/2
        extra_args = tuple([center_idx])
        
    else:
        sys.exit('Could not find filtering function for alias: %s' % function)
    
    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    shutil.copy2(params, out_dir)
        
    print '\nReading input raster...\n'
    t1 = time.time()
    ds = gdal.Open(path)
    band = ds.GetRasterBand(databand)
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    driver = ds.GetDriver()
    
    # Get an array and mask out nodata values with nans
    if 'nodata' not in inputs:
        print 'nodata not specified in params. Getting nodata value from input dataset...\n'
        nodata = band.GetNoDataValue()
    ar = band.ReadAsArray()
    ds = None
    array_dtype = ar.dtype
    ar = ar.astype(np.float32)
    mask = (ar != nodata) & (ar != 255)
    ar[~mask] = np.nan
    if 'area' in function.lower():
        ar[(ar != filter_value) & mask] = 0
    #import pdb; pdb.set_trace()
    ysize, xsize = ar.shape
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    if kernel_type.lower() == 'circle':
        #kernel_size /= 2
        kernel = circle_mask(kernel_size)
    else:
        kernel = np.ones((kernel_size, kernel_size))
    tile_buffer = kernel.shape[0]/2
    # Tile up the array to filter in parallel
    # Find empty tiles
    print 'Finding empty tiles...'
    t1 = time.time()
    df_tiles, df_tiles_rc, _ = get_tiles(n_tiles, xsize, ysize, tx)

    total_tiles = len(df_tiles)
    empty_tiles = find_empty_tiles(df_tiles, mask, tx)
    df_tiles = df_tiles_rc.select(lambda x: x not in empty_tiles)
    print '%s empty tiles found of %s total tiles\n%.1f minutes\n' %\
    (len(empty_tiles), total_tiles, (time.time() - t1)/60)
    
    # Add buffer around each tile
    df_buf = df_tiles.copy()
    df_buf[['ul_r', 'ul_c']] = df_buf[['ul_r', 'ul_c']] - tile_buffer
    df_buf[['lr_r', 'lr_c']] = df_buf[['lr_r', 'lr_c']] + tile_buffer
    df_buf[['ul_r', 'lr_r']] = df_buf[['ul_r', 'lr_r']].clip(0, ysize)
    df_buf[['ul_c', 'lr_c']] = df_buf[['ul_c', 'lr_c']].clip(0, xsize)
    
    # Get arrays
    print 'Getting buffered arrays...'
    t1 = time.time()
    n_full_tiles = len(df_tiles)
    args = []
    for i, (ind, r) in enumerate(df_buf.iterrows()):
        this_ar = ar[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        args.append([ind, this_ar, func, kernel, extra_args, i + 1, n_full_tiles])
        #arrays.append([i, this_ar])
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    print 'Filtering chunks in parallel with %s jobs...' % n_jobs
    p = Pool(n_jobs)
    tiles = p.map(par_filter, args, 1)
    print '\nTotal time for filtering: %.1f minutes\n' % ((time.time() - t1)/60)
    #tiles = arrays'''
    
    print 'Tiling pieces back together...'
    t1 = time.time()
    filtered = np.full(ar.shape, nodata, dtype=array_dtype)
    #tiles = [[i,np.ones((df_buf.ix[60,'lr_r'] - df_buf.ix[60,'ul_r'], df_buf.ix[60,'lr_c'] - df_buf.ix[60,'ul_c']))]]
    for i, buffered_tile in tiles:
        b_inds = df_buf.ix[i, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        t_inds = df_tiles.ix[i, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        d_ulr, d_lrr, d_ulc, d_lrc = t_inds - b_inds
        
        tile = buffered_tile[d_ulr : d_lrr, d_ulc : d_lrc]
        tile[np.isnan(tile)] = nodata
        tile = tile.astype(array_dtype)
        t_ulr, t_lrr, t_ulc, t_lrc = t_inds
        filtered[t_ulr : t_lrr, t_ulc : t_lrc] = tile
    print '%.1f minutes\n' % ((time.time() - t1)/60)   
    
    #filtered = filtered.astype(array_dtype)
    if 'out_nodata' in inputs: nodata = int(inputs['out_nodata'])
    filtered[np.isnan(filtered) | ~mask] = nodata
    
    if ar.max() <= 255 and ar.min() >= 0:
        gdal_dtype = gdal.GDT_Byte
    else:
        gdal_dtype = gdal.GDT_UInt16
    array_to_raster(filtered, tx, prj, driver, out_path, gdal_dtype, nodata)
    del ar, filtered, tiles, args, p
    
    print 'Total time: %.1f minutes' % ((time.time() - t0)/60)


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))
    
    
    