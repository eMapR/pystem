# -*- coding: utf-8 -*-
"""
Align one raster to the exact extent of another, filling the input raster with 
nodata values if necessary. If a mask value is given, the input raster will also
be masked where snap_raster == mask_val.


"""
import os
import sys
import gdal
import time
import shutil
import numpy as np
import cPickle as pickle

package_dir = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(package_dir)
import stem
from mosaic_by_tsa import get_offset_array_indices, calc_offset, get_min_numpy_dtype
from lthacks import createMetadata, array_to_raster, get_gdal_driver


def make_tiles(n_tiles, ds_snap):
    
    try:
        n_tiles = int(n_tiles)
    except ValueError:
        if ',' in n_tiles:
            try: n_tiles = [int(i) for i in n_tiles.split(',')]
            except ValueError: pass
        else:
            try: n_tiles =[int(i) for i in n_tiles.split()]
            except ValueError:
                raise ValueError('format of n_tiles not understood: %s' % n_tiles)
        
    ysize = ds_snap.RasterYSize
    xsize = ds_snap.RasterXSize
    tx = ds_snap.GetGeoTransform()
    
    # Figure out how many tiles belong in each row if necessary
    if isinstance(n_tiles, int):
        # if n_tiles = nx * ny and nx = ny * ratio -> y = (n_tiles/ratio) ** .5
        ratio = xsize/float(ysize)
        ny = int((n_tiles/ratio) ** .5)
        nx = int(n_tiles/ny)
        n_tiles = ny, nx
    
    _, tiles, __ = stem.get_tiles(n_tiles, xsize, ysize, tx)
    #stem.coords_to_shp(_, '/vol/v2/stem/extent_shp/CAORWA.shp', '/home/server/pi/homes/shooper/delete/tiles.shp')
    
    return tiles
        

def snap_array(ds_in, ds_snap, tx_in, tx_snap, in_nodata, out_nodata, mask_val=None):
    
    ar_in = ds_in.ReadAsArray()
    if mask_val is not None:
        ar_snap = ds_snap.ReadAsArray()
    in_shape = ar_in.shape
    out_shape = ds_snap.RasterYSize, ds_snap.RasterXSize
    
    offset = calc_offset((tx_snap[0], tx_snap[3]), tx_in)
    snap_inds, in_inds = get_offset_array_indices(out_shape, in_shape, offset)
    np_dtype = ar_in.dtype
    ar = np.full(out_shape, out_nodata, dtype=np_dtype)
    ar_in[ar_in == in_nodata] = out_nodata
    ar[snap_inds[0]:snap_inds[1], snap_inds[2]:snap_inds[3]] = ar_in[in_inds[0]:in_inds[1], in_inds[2]:in_inds[3]]
    
    if mask_val is not None:
        mask_val = int(mask_val)
        ar[ar_snap == mask_val] = out_nodata
    
    return ar


def snap_by_tile(ds_in, ds_snap, tiles, tx_snap, tx_in, in_nodata, out_nodata, out_dir, mask_val=None):
    
    prj = ds_in.GetProjection()
    driver = gdal.GetDriverByName('gtiff')
    
    if mask_val is not None:
        mask_val = int(mask_val)
        
    row_off, col_off = calc_offset((tx_snap[0], tx_snap[3]), tx_in)
    in_size = ds_in.RasterYSize, ds_in.RasterXSize
    
    n_tiles = float(len(tiles))
    t1 = time.time()
    msg = '\rProccessing tile %d/%d (%.1f%%) || %.1f/~%.1f minutes'
    
    template = os.path.join(out_dir, 'tile_%s.pkl')
    mins = []
    maxs = []
    for i, (tile_id, coords) in enumerate(tiles.iterrows()):
        
        tile_off = row_off - coords.ul_r, col_off - coords.ul_c
        tile_size = coords.lr_r - coords.ul_r, coords.lr_c - coords.ul_c
        tile_inds, in_inds = get_offset_array_indices(tile_size, in_size, tile_off)
        
        in_ulr, in_lrr, in_ulc, in_lrc = in_inds
        in_xsize = in_lrc - in_ulc
        in_ysize = in_lrr - in_ulr
        if in_xsize <= 0 or in_ysize <= 0: # They don't overlap
            continue
        ar_in = ds_in.ReadAsArray(in_ulc, in_ulr, in_xsize, in_ysize)
        if np.all(ar_in == in_nodata):
            continue
        ar_out = np.full(tile_size, out_nodata, dtype=ar_in.dtype)
        ar_out[tile_inds[0]:tile_inds[1], tile_inds[2]:tile_inds[3]] = ar_in
        ar_out[ar_out == in_nodata] = out_nodata
        if mask_val is not None:
            mask = ds_snap.ReadAsArray(coords.ul_c, coords.ul_r, tile_size[1], tile_size[0]) == mask_val
            ar_out[mask] = out_nodata
        
        out_path = template % tile_id
        with open(out_path, 'wb') as f:
            pickle.dump(ar_out, f, protocol=-1)
        mins.append(ar_out.min())
        maxs.append(ar_out.max())
        tiles.loc[tile_id, 'file'] = out_path
        
        cum_time = (time.time() - t1)/60.
        est_time = cum_time/(i + 1) * (n_tiles - i) # estimate remaing time
        sys.stdout.write(msg % (i + 1, n_tiles, (i + 1)/n_tiles * 100, cum_time, est_time))
        sys.stdout.flush()
        
        
        '''ulx, xres, _, uly, _, yres = tx_snap
        tx = coords.ul_c * xres + ulx, xres, 0, coords.ul_r * yres + uly, 0, yres
        array_to_raster(ar_out, tx, prj, driver, '/home/server/pi/homes/shooper/delete/tile_%s.tif' % tile_id, gdal.GDT_Int16, out_nodata)'''
        
    dtype = get_min_numpy_dtype(np.array(mins + maxs))
    
    return dtype
    

def main(in_raster, snap_raster, in_nodata, out_nodata, out_path=None, mask_val=None, overwrite=False, n_tiles=None):
    
    t0 = time.time()
    in_nodata = int(in_nodata)
    out_nodata = int(out_nodata)
    
    print '\nOpening datasets... '
    t1 = time.time()
    ds_in = gdal.Open(in_raster)
    tx_in = ds_in.GetGeoTransform()

    ds_snap = gdal.Open(snap_raster)
    snap_size = ds_snap.RasterYSize, ds_snap.RasterXSize
    tx_snap = ds_snap.GetGeoTransform()
    prj = ds_snap.GetProjection()
    print '%.1f seconds\n' % (time.time() - t1)
    
    if n_tiles:
        if not out_path:
            raise IOError('n_tiles was given, but no out_path specified')
        tiles = make_tiles(n_tiles, ds_snap)
        temp_dir = os.path.join(os.path.dirname(out_path), 'temp_tiles')
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
        np_dtype = snap_by_tile(ds_in, ds_snap, tiles, tx_snap, tx_in, in_nodata, out_nodata, temp_dir, mask_val)
        ar = np.full(snap_size, out_nodata, dtype=np_dtype)
        for tile_id, coords in tiles.dropna(subset=['file']).iterrows():
            with open(coords.file, 'rb') as f:
                ar_tile = pickle.load(f)
            t_ysize, t_xsize = ar_tile.shape
            ul_r = coords.ul_r
            ul_c = coords.ul_c
            lr_r = ul_r + t_ysize
            lr_c = ul_c + t_xsize
            ar[ul_r : lr_r, ul_c : lr_c] = ar_tile
        shutil.rmtree(temp_dir)
        
    else:
        ar_in = ds_in.ReadAsArray()
        ar_snap = ds_snap.ReasAsArray()
        ar = snap_array(ar_in, ar_snap, tx_in, tx_snap, in_nodata, out_nodata, mask_val)
        print '%.1f minutes\n' % ((time.time() - t1)/60)
    
    if out_path:
        if ar.max() <= 255 and ar.min() >= 0:
            gdal_dtype = gdal.GDT_Byte
        else:
            gdal_dtype = gdal.GDT_Int16
        
        if os.path.exists(out_path) and not overwrite: 
            sys.exit('out_path already exists')
        driver = get_gdal_driver(out_path)
        array_to_raster(ar, tx_snap, prj, driver, out_path, gdal_dtype, out_nodata)
        
        # Write metadata
        desc = ('Input raster %s snapped to the extent of %s.') % (in_raster, snap_raster)
        if mask_val:
            desc += ' Data were masked from snap raster with value %s.' % mask_val
        createMetadata(sys.argv, out_path, description=desc)
    else:
        return ar
    
    print '\nTotal time to snap raster: %.1f seconds\n' % (time.time() - t0)


if __name__ == '__main__':
    
    sys.exit(main(*sys.argv[1:]))
    
