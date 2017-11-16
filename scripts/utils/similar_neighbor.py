# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 09:20:02 2017

@author: shooper
"""
import os
import sys
import time
import numpy as np
from osgeo import gdal
from multiprocessing import Pool

from evaluation import get_samples
from stem import get_tiles, find_empty_tiles
from lthacks import createMetadata
from lthacks import array_to_raster

def buffered_tile_inds(n_tiles, xsize, ysize, tx, tile_buffer, mask):
    
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
    
    return df_tiles, df_buf


def par_get_match(args):
    
    t0 = time.time()
    tile_ind, this_in, this_match, in_nodata, match_nodata, count, total_tiles = args
    matched_vals, _ = get_samples(this_in, this_match, in_nodata, match_nodata, match=True)
    print 'Time for getting array %s of %s: %.1f seconds' % (count, total_tiles, time.time() - t0)
    return tile_ind, matched_vals
    

def main(in_raster, match_raster, in_nodata, match_nodata, out_raster=None, n_jobs=25):
    
    t0 = time.time()
    in_nodata = int(in_nodata)
    match_nodata = int(match_nodata)
    
    print '\nReading datasets... '
    t1 = time.time()
    ds_in = gdal.Open(in_raster)
    ar_in = ds_in.ReadAsArray()
    tx = ds_in.GetGeoTransform()
    prj = ds_in.GetProjection()
    driver = ds_in.GetDriver()
    ds_in = None
    
    ds_match = gdal.Open(match_raster)
    ar_match = ds_match.ReadAsArray()
    ds_match = None
    print '%.1f seconds\n' % (time.time() - t1)
    
    n_tiles = (30, 20)
    ysize, xsize = ar_in.shape
    tile_buffer = 1
    nodata_mask = (ar_in != in_nodata) & (ar_match != match_nodata)
    tile_inds, buf_inds = buffered_tile_inds(n_tiles, xsize, ysize, tx, tile_buffer, nodata_mask)
    total_tiles = len(tile_inds)
    
    args = []
    for i, (ind, r) in enumerate(buf_inds.iterrows()):
        this_in = ar_in[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        this_match = ar_match[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        args.append([ind, this_in, this_match, in_nodata, match_nodata, i+1, total_tiles])
    
    #n_jobs = 10
    print 'Filtering chunks in parallel with %s jobs...' % n_jobs
    p = Pool(n_jobs)
    matched_arrays = p.map(par_get_match, args, 1)
    print '\nTotal time for filtering: %.1f minutes\n' % ((time.time() - t1)/60)
    
    ar = np.full((ysize, xsize), in_nodata, dtype=ar_in.dtype)
    for ind, matched_vals in matched_arrays:
        r = buf_inds.ix[ind]
        this_ysize = r.lr_r - r.ul_r
        this_xsize = r.lr_c - r.ul_c
        buffered_tile = np.full((this_ysize, this_xsize), in_nodata, dtype=ar_in.dtype)
        this_mask = nodata_mask[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        buffered_tile[this_mask] = matched_vals
        
        b_inds = buf_inds.ix[ind, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        t_inds = tile_inds.ix[ind, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        d_ulr, d_lrr, d_ulc, d_lrc = t_inds - b_inds
        
        tile = buffered_tile[d_ulr : d_lrr, d_ulc : d_lrc]
        t_ulr, t_lrr, t_ulc, t_lrc = t_inds
        ar[t_ulr : t_lrr, t_ulc : t_lrc] = tile
        
    '''ar = np.full((ysize, xsize), in_nodata, dtype=ar_in.dtype)
    for i, (ind, r) in enumerate(buf_inds.iterrows()):
        #in_arrays.append([ind, ar_in[r.ul_r : r.lr_r, r.ul_c : r.lr_c]])
        #match_arrays.append([ind, ar_match[r.ul_r : r.lr_r, r.ul_c : r.lr_c]])
        print 'Getting array %s of %s...' % (i + 1, total_tiles)
        t2 = time.time()
        this_in = ar_in[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        this_match = ar_match[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        this_mask = nodata_mask[r.ul_r : r.lr_r, r.ul_c : r.lr_c]
        
        _, matched_vals = get_samples(this_match, this_in, in_nodata, match_nodata, match=True)
        this_ysize = r.lr_r - r.ul_r
        this_xsize = r.lr_c - r.ul_c
        buffered_tile = np.full((this_ysize, this_xsize), in_nodata, dtype=ar_in.dtype)
        buffered_tile[this_mask] = matched_vals
        
        b_inds = buf_inds.ix[ind, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        t_inds = tile_inds.ix[ind, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        d_ulr, d_lrr, d_ulc, d_lrc = t_inds - b_inds
        
        tile = buffered_tile[d_ulr : d_lrr, d_ulc : d_lrc]
        #tile[np.isnan(tile)] = nodata
        #tile = tile.astype(array_dtype)
        t_ulr, t_lrr, t_ulc, t_lrc = t_inds
        #import pdb; pdb.set_trace()
        ar[t_ulr : t_lrr, t_ulc : t_lrc] = tile
        print '%.1f seconds\n' % (time.time() - t2)'''
    
    '''t1 = time.time()
    ar = np.full((xsize, ysize), in_nodata, dtype=np.int16)
    for i, ind in enumerate(buf_inds.index):
        t2 = time.time()
        print 'Getting array %s of %s...'
        matched_vals, _ = get_samples(ar_in, ar_match, in_nodata, match_nodata, match=True)
        buffered_tile = 
        b_inds = buf_inds.ix[i, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        t_inds = tile_inds.ix[i, ['ul_r', 'lr_r', 'ul_c', 'lr_c']]
        d_ulr, d_lrr, d_ulc, d_lrc = t_inds - b_inds
        
        tile = buffered_tile[d_ulr : d_lrr, d_ulc : d_lrc]
        tile[np.isnan(tile)] = nodata
        tile = tile.astype(array_dtype)
        t_ulr, t_lrr, t_ulc, t_lrc = t_inds
        filtered[t_ulr : t_lrr, t_ulc : t_lrc] = tile
        print '%.1f seconds\n' % (time.time() - t2)
    
    matched_vals, _ = get_samples(ar_in, ar_match, in_nodata, match_nodata, match=True)
    ar = np.full((xsize, ysize), in_nodata, dtype=np.int16)
    ar[nodata_mask] = matched_vals
    print '%.1f seconds\n' % (time.time() - t1)
    '''    
    if not out_raster:
        in_dir, in_bn = os.path.split(in_raster)
        in_name, ext = os.path.splitext(in_bn)
        match_bn, _ = os.path.splitext(os.path.basename(match_raster))
        out_bn = '{0}_matched_{1}{2}'.format(in_name, match_bn, ext)
        out_raster = os.path.join(in_dir, out_bn)

    if ar.max() <= 255 and ar.min() >= 0:
        gdal_dtype = gdal.GDT_Byte
    else:
        gdal_dtype = gdal.GDT_Int16
    array_to_raster(ar, tx, prj, driver, out_raster, gdal_dtype, match_nodata)
    
    # Write metadata
    desc = ('Raster with values where input raster %s matched another ' +\
            'raster %s best within a 3x3 kernel.') % (in_raster, match_raster)
    createMetadata(sys.argv, out_raster, description=desc)
    
    print '\nTotal time to match raster: %.1f minutes\n' % ((time.time() - t0)/60)
    

if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))