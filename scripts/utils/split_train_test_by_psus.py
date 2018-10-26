# -*- coding: utf-8 -*-
"""
Created on Fri May 11 12:08:24 2018

@author: shooper

eample command:
python split_train_test_by_psus.py  /vol/v1/general_files/datasets/spatial_data/nlcd/nlcd_2001_v2/nlcd_2001_landcover_clipped_to_conus_tiles.tif 255 /vol/v2/stem/conus/vector/sampling_psu_100_2001.shp /vol/v2/stem/conus/reference_rasters
"""

import os, sys
from osgeo import gdal, gdal_array
import pandas as pd
import numpy as np

package_dir = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(package_dir)
from stem import calc_offset

from lthacks import attributes_to_df, array_to_raster, createMetadata


def main(raster, nodata, psu_shp, out_dir):
    
    nodata = int(nodata)
    psus = attributes_to_df(psu_shp)
    ds = gdal.Open(raster)
    ar = ds.GetVirtualMemArray()#ReadAsArray()
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    driver = gdal.GetDriverByName('gtiff')
    
    # Just extract the test sample first
    test_sample_dfs = []
    print '\nGetting test samples for PSUs...'
    for i, psu in psus.iterrows():
        # Calc offsets
        row_off, col_off = calc_offset((tx[0], tx[3]), psu[['ul_x', 'ul_y']], tx)
        n_rows = abs(int((psu.ymax - psu.ymin)/tx[5]))
        n_cols = abs(int((psu.xmax - psu.xmin)/tx[1]))
        
        # Get values
        test_data = ar[row_off : row_off + n_rows, col_off : col_off + n_cols].ravel()
        mask = test_data != nodata
        
        # Get row/col and x/y vals
        test_data = test_data[mask]
        row_inds, col_inds = np.indices((n_rows, n_cols), dtype=np.uint32)
        row_inds = row_inds.ravel()[mask]
        col_inds = col_inds.ravel()[mask]
        these_row_inds = row_inds + row_off
        these_col_inds = col_inds + col_off
        y_coords = row_inds * tx[5] + psu.ul_y
        x_coords = col_inds * tx[1] + psu.ul_x

        df = pd.DataFrame({'row': these_row_inds,
                           'col': these_col_inds,
                           'y': y_coords,
                           'x': x_coords,
                           'value': test_data,
                           'tile_id': psu['name']
                           })
        #import pdb; pdb.set_trace()
        test_sample_dfs.append(df)
        
    test_sample = pd.concat(test_sample_dfs, ignore_index=True)
    basename = os.path.basename(raster)
    out_txt = os.path.join(out_dir, basename.replace(basename[-4:], '_test.txt'))
    test_sample.to_csv(out_txt, sep='\t', index=False)
    
    # Read the raster as a write-able array and set all test samples to nodata
    print '\nAssigning nodata val to PSUs in training raster...\n'
    ar = ds.ReadAsArray()
    ar[test_sample.row, test_sample.col] = nodata
    out_raster = out_txt.replace('_test.txt', '_train.tif')
    array_to_raster(ar, tx, prj, driver, out_raster, nodata=nodata)
    
    desc = 'Training raster and test sample (text file with the same name but "_test" at the end) for making and evaluating STEM CONUS maps. Primary sampling units (PSUs) reserved for testing are assigned nodata.'
    desc += '\n\tInput raster: %s' % os.path.abspath(raster)
    desc += '\n\tNodata value: %s' % nodata
    desc += '\n\tPSU shapefile: %s' % os.path.abspath(psu_shp)
    desc += '\n\tOutput directory: %s\n' % os.path.abspath(out_dir)
    createMetadata(sys.argv, out_raster, description=desc)#'''
    
    ds = None
    
if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
        