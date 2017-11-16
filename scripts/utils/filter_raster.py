# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:13:57 2016

@author: shooper
"""
import os
import sys
import time
import shutil
import numpy as np
from osgeo import gdal, gdalnumeric
from multiprocessing import Pool
from scipy import ndimage as ndi

from stem import get_tiles, find_empty_tiles, mode, coords_to_shp
from mosaic_by_tsa import array_to_raster

from lthacks import createMetadata, write_params_to_meta

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
    #ind, ar, func, kernel, extra_args, i, n_tiles = args
    ind, path, databand, nodata, r, tile_coords, out_dir, func, kernel, extra_args, i, n_tiles = args
    ds = gdal.Open(path)
    nrows = r.lr_r - r.ul_r
    ncols = r.lr_c - r.ul_c
    ar = ds.GetRasterBand(databand).ReadAsArray(r.ul_c, r.ul_r, ncols, nrows)
    mask = ar == nodata
    if np.all(mask):
        return ind, None
    ar = ndi.generic_filter(ar, func, footprint=kernel, extra_arguments=extra_args)
    ar[mask] = nodata
    
    _, x_res, _, _, _, y_res = ds.GetGeoTransform()
    driver = gdal.GetDriverByName('gtiff')
    prj = ds.GetProjection()
    tx = tile_coords.ul_x, x_res, 0, tile_coords.ul_y, 0, y_res
    out_path = os.path.join(out_dir, 'tile_%s.tif' % ind)
    array_to_raster(ar, tx, prj, driver, out_path, nodata=nodata)
    
    print 'Time for tile %s of %s: %.1f minutes' % (i, n_tiles, ((time.time() - t0)/60))
    ds = None
    
    return ind, out_path


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
    if 'kernel_type' in inputs: kernel_type = inputs['kernel_type']
    
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
    # shutil.copy2(params, out_dir) # don't need to copy because params are written to meta
        
    print '\nReading input raster...\n'
    t1 = time.time()
    ds = gdal.Open(path)
    band = ds.GetRasterBand(databand)
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    driver = ds.GetDriver()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    
    # Get an array and mask out nodata values with nans
    if 'nodata' not in inputs:
        print 'nodata not specified in params. Getting nodata value from input dataset...\n'
        nodata = band.GetNoDataValue()
    '''ar = band.ReadAsArray()
    ds = None
    array_dtype = ar.dtype
    ar = ar.astype(np.float16)
    mask = (ar != nodata) #& (ar != 255)
    ar[~mask] = np.nan'''
    if 'area' in function.lower():
        ar[(ar != filter_value) & mask] = 0
    #import pdb; pdb.set_trace()
    #ysize, xsize = ar.shape
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
    '''empty_tiles = find_empty_tiles(df_tiles, mask, tx)
    df_tiles = df_tiles_rc.select(lambda x: x not in empty_tiles)
    print '%s empty tiles found of %s total tiles\n%.1f minutes\n' %\
    (len(empty_tiles), total_tiles, (time.time() - t1)/60)'''
    
    # Add buffer around each tile
    df_buf = df_tiles_rc.copy()
    df_buf[['ul_r', 'ul_c']] = df_buf[['ul_r', 'ul_c']] - tile_buffer
    df_buf[['lr_r', 'lr_c']] = df_buf[['lr_r', 'lr_c']] + tile_buffer
    df_buf[['ul_r', 'lr_r']] = df_buf[['ul_r', 'lr_r']].clip(0, ysize)
    df_buf[['ul_c', 'lr_c']] = df_buf[['ul_c', 'lr_c']].clip(0, xsize)
    
    # Get arrays
    print 'Getting buffered arrays...'
    t1 = time.time()
    n_full_tiles = len(df_tiles)
    args = []
    temp_dir = os.path.join(out_dir, 'tiles')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    for i, (ind, r) in enumerate(df_buf.iterrows()):
        #this_ar = ar[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        #args.append([ind, this_ar, func, kernel, extra_args, i + 1, n_full_tiles])
        args.append([ind, path, databand, nodata, r, df_tiles.ix[ind], temp_dir, func, kernel, extra_args, i + 1, n_full_tiles])
        #arrays.append([i, this_ar])
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    print 'Filtering chunks in parallel with %s jobs...' % n_jobs
    p = Pool(n_jobs)
    tiles = p.map(par_filter, args, 1)

    print '\nTotal time for filtering: %.1f minutes\n' % ((time.time() - t1)/60)#'''

    
    print 'Tiling pieces back together...'
    t1 = time.time()
    gdal_dtype = band.DataType
    array_dtype = gdalnumeric.GDALTypeCodeToNumericTypeCode(gdal_dtype)
    filtered = np.full((ysize, xsize), nodata, dtype=array_dtype)
    for i, tile_path in tiles:
        if not tile_path:
            continue
        ds_t = gdal.Open(tile_path)
        buffered_tile = ds_t.ReadAsArray()
        b_inds = df_buf.ix[i, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        t_inds = df_tiles_rc.ix[i, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        d_ulr, d_lrr, d_ulc, d_lrc = t_inds - b_inds
        tile = buffered_tile[d_ulr : d_lrr, d_ulc : d_lrc]
        tile[np.isnan(tile)] = nodata
        tile = tile.astype(array_dtype)
        t_ulr, t_lrr, t_ulc, t_lrc = t_inds
        filtered[t_ulr : t_lrr, t_ulc : t_lrc] = tile
    print '%.1f minutes\n' % ((time.time() - t1)/60)   
    
    #filtered = filtered.astype(array_dtype)
    if 'out_nodata' in inputs: 
        #filtered[np.isnan(filtered) | ~mask] = nodata
        filtered[filtered == nodata] = int(inputs['out_nodata'])
        nodata = int(inputs['out_nodata'])

    try:
        array_to_raster(filtered, tx, prj, driver, out_path, dtype=gdal_dtype, nodata=nodata)
    except:
        array_to_raster(filtered, tx, prj, driver, out_path, gdal.GDT_Byte, nodata=nodata)
    desc = ('Raster filtered by kernel of shape {kernel_type} and size ' +\
            '{kernel_size} and function {func}').format(kernel_type=kernel_type,
                                                        kernel_size=kernel_size, 
                                                        func=function)
    meta_path = createMetadata(sys.argv, out_path, description=desc)
    write_params_to_meta(meta_path, params)
    del ar, filtered, tiles, args, p
    ds = None
    import pdb; pdb.set_trace()
    shutil.rmtree(temp_dir)
    
    print 'Total time: %.1f minutes' % ((time.time() - t0)/60)


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))
    
    
    