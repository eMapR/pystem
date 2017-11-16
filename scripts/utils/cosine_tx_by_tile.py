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


def transform_aspect(aspect_raster, nodata, transform='both', out_dir=None):
    
    nesw = False
    nwse = False
    
    ds = gdal.Open(aspect_raster)
    ar = ds.ReadAsArray()
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    ds = None
    
    driver = gdal.GetDriverByName('envi')
    if not out_dir:
        out_dir = os.path.dirname(aspect_raster)
    
    if transform == 'both':
        nesw = True
        nwse = True
    elif transform == 'nesw':
        nesw = True
    else:
        nwse = True
    mask = ar != nodata
    if nesw:
        ar_nesw = np.full(ar.shape, 255, dtype=np.uint8)
        ar_nesw[mask] = 100 * (np.cos(np.radians(ar - 225)) + 1)[mask]
        out_path = aspect_raster.replace('.bsq', '_nesw.bsq') 
        array_to_raster(ar_nesw, tx, prj, driver, out_path, dtype=gdal.GDT_Byte, silent=True)
        desc = ('Cosine transformed aspect representing northeast/southwest-ness.\n' +\
            '\tCreated with the expression: 100 * (np.cos(np.radians(ar - 225)) + 1)\n' + \
            '\tInput aspect raster: {0}').format(aspect_raster)
        createMetadata(sys.argv, out_path, description=desc)
    if nwse:
        ar_nwse = np.full(ar.shape, 255, dtype=np.uint8)
        ar_nwse[mask] = 100 * (np.cos(np.radians(ar - 135)) + 1)[mask]
        out_path = aspect_raster.replace('.bsq', '_nwse.bsq') 
        array_to_raster(ar_nwse, tx, prj, driver, out_path, dtype=gdal.GDT_Byte, silent=True)
        desc = ('Cosine transformed aspect representing northwest/southeast-ness.\n' +\
            '\tCreated with the expression: 100 * (np.cos(np.radians(ar - 135)) + 1)\n' + \
            '\tInput aspect raster: {0}').format(aspect_raster)
        createMetadata(sys.argv, out_path, description=desc)


def par_transform_aspect(args):
    
    t0 = time.time()
    (aspect_raster, nodata, transform, out_dir, n), n_tiles = args
    transform_aspect(aspect_raster, nodata, transform, out_dir)
    print 'Transformed aspect for %s of %s tiles: %.1f seconds' % (n, n_tiles, (time.time() -t0))


def main(search_dir, search_str, nodata, out_dir=None, fileList=None, transform='both'):
    
    t0 = time.time()
    nodata = int(nodata)
    
    i = 1
    args = []
    for root, dirs, files in os.walk(search_dir):
        for f in fnmatch.filter(files, search_str):
            path = os.path.join(root, f)
            args.append([path, nodata, transform, out_dir, i])
            i += 1
    n_tiles = len(args)
    args = [[a, n_tiles] for a in args]
    
    pool = Pool(20)
    pool.map(par_transform_aspect, args, 1)#'''
    
    print '\nFinished in %.1f minutes\n' % ((time.time() - t0)/60)


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
    
    
            