# -*- coding: utf-8 -*-
"""
Created on Mon May  1 10:04:05 2017

@author: shooper
"""

import os
import sys
import time
import fnmatch
import numpy as np
from multiprocessing import Pool
from osgeo import gdal

from lthacks import array_to_raster, createMetadata


def scale_raster(in_raster, nodata, scale):
    
    ds = gdal.Open(in_raster)
    ar = ds.ReadAsArray()
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    driver = ds.GetDriver()
    ds = None
    
    #driver = gdal.GetDriverByName('envi')
    mask = ar != nodata
    ar[mask] = ar[mask] * scale
    out_path = in_raster#.replace('.bsq', '_scaled%s.bsq' % scale) 
    array_to_raster(ar, tx, prj, driver, out_path, dtype=gdal.GDT_Int16, silent=True)
    desc = '%s scaled by %s.\n' % (in_raster, scale)
    createMetadata(sys.argv, out_path, description=desc)


def par_scale_raster(args):
    
    t0 = time.time()
    (in_raster, nodata, scale, n), n_tiles = args
    scale_raster(in_raster, nodata, scale)
    print 'Transformed aspect for %s of %s tiles: %.1f seconds' % (n, n_tiles, (time.time() -t0))


def main(search_dir, search_str, nodata, scale, out_dir=None):
    
    t0 = time.time()
    nodata = int(nodata)
    scale = int(scale)
    
    i = 1
    args = []
    for root, dirs, files in os.walk(search_dir):
        for f in fnmatch.filter(files, search_str):
            in_raster = os.path.join(root, f)
            args.append([in_raster, nodata, scale, i])
            i += 1
    n_tiles = len(args)
    args = [[a, n_tiles] for a in args]
    
    pool = Pool(20)
    pool.map(par_scale_raster, args, 1)#'''
    
    print '\nFinished in %.1f minutes\n' % ((time.time() - t0)/60)


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
    
    
            