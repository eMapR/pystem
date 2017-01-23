# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 10:10:19 2016

@author: shooper
"""

import sys
import gdal
import time
import numpy as np
from mosaic_by_tsa import get_offset_array_indices, calc_offset, array_to_raster


def main(in_raster, snap_raster, in_nodata, out_nodata, out_path=None):
    
    t0 = time.time()
    in_nodata = int(in_nodata)
    out_nodata = int(out_nodata)
    
    print '\nOpening datasets... '
    t1 = time.time()
    ds_in = gdal.Open(in_raster)
    ar_in = ds_in.ReadAsArray()
    tx_in = ds_in.GetGeoTransform()
    driver = ds_in.GetDriver()
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
    in_nodata = ar_in.min()
    ar_in[ar_in == in_nodata] = out_nodata
    ar[snap_inds[0]:snap_inds[1], snap_inds[2]:snap_inds[3]] = ar_in[in_inds[0]:in_inds[1], in_inds[2]:in_inds[3]]
    ar[ar > 60000] = out_nodata - 1
    print '%.1f seconds\n' % (time.time() - t1)
    
    if ar.max() <= 255 and ar.min() >= 0:
        gdal_dtype = gdal.GDT_Byte
    else:
        gdal_dtype = gdal.GDT_UInt16
    
    if not out_path: out_path = in_raster
    array_to_raster(ar * 10, tx_snap, prj, driver, out_path, gdal_dtype, out_nodata)
    
    print '\nTotal time to snap raster: %.1f seconds\n' % (time.time() - t0)


if __name__ == '__main__':
    
    sys.exit(main(*sys.argv[1:]))