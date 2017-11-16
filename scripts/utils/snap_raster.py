# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 10:10:19 2016

@author: shooper
"""
import os
import sys
import gdal
import time
import numpy as np

from mosaic_by_tsa import get_offset_array_indices, calc_offset
from lthacks import createMetadata, array_to_raster, get_gdal_driver


def main(in_raster, snap_raster, in_nodata, out_nodata, out_path=None, mask_val=None, overwrite=False):
    
    t0 = time.time()
    in_nodata = int(in_nodata)
    out_nodata = int(out_nodata)
    
    print '\nOpening datasets... '
    t1 = time.time()
    ds_in = gdal.Open(in_raster)
    ar_in = ds_in.ReadAsArray()
    tx_in = ds_in.GetGeoTransform()
    #driver = ds_in.GetDriver()
    ds_in = None
    
    ds_snap = gdal.Open(snap_raster)
    ar_snap = ds_snap.ReadAsArray()
    tx_snap = ds_snap.GetGeoTransform()
    prj = ds_snap.GetProjection()
    ds_snap = None
    print '%.1f seconds\n' % (time.time() - t1)
    
    print 'Snapping input raster...'
    t1 = time.time()
    offset = calc_offset((tx_snap[0], tx_snap[3]), tx_in)
    snap_inds, in_inds = get_offset_array_indices(ar_snap.shape, ar_in.shape, offset)
    np_dtype = ar_in.dtype
    ar = np.full(ar_snap.shape, out_nodata, dtype=np_dtype)
    ar_in[ar_in == in_nodata] = out_nodata
    ar[snap_inds[0]:snap_inds[1], snap_inds[2]:snap_inds[3]] = ar_in[in_inds[0]:in_inds[1], in_inds[2]:in_inds[3]]
    
    if mask_val:
        mask_val = int(mask_val)
        ar[ar_snap == mask_val] = out_nodata
    
    print '%.1f seconds\n' % (time.time() - t1)
    
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
    
