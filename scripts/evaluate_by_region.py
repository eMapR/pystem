# -*- coding: utf-8 -*-
"""
Created on Thu May 11 11:53:45 2017

@author: shooper
"""

import os
import sys
import time
import random
import shutil
from osgeo import gdal
import pandas as pd
import numpy as np

import evaluation
import stem
from get_stratified_random_pixels import parse_bins, read_params

from lthacks import stats_functions as sf
from lthacks import createMetadata, df_to_shp

'''n_samples
train_txt
region_ids
bins
p_nodata
t_nodata
p_raster
t_raster
region_raster
out_dir'''

def main(params, test_txt=None, region_nodata=0, region_shp=None):
    
    t0 = time.time()
    # Read params and make variables from each line
    inputs = read_params(params)
    for var in inputs:
        exec("{0} = str({1})").format(var, inputs[var])
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    shutil.copy2(params, out_dir)
    
    n_samples = int(n_samples)
    p_nodata = int(p_nodata)
    t_nodata = int(t_nodata)
    bins = parse_bins(bins)
    if 'region_nodata' in inputs:
        region_nodata = int(region_nodata)
    
    ds = gdal.Open(region_raster)
    regions = ds.ReadAsArray()
    nodata_mask = regions == region_nodata
    ds = None
    
    ds_t = gdal.Open(t_raster)
    #ar_t = ds_t.ReadAsArray()
    
    ds_p = gdal.Open(p_raster)
    #ar_p = ds_p.ReadAsArray()
    d = {}
    if 'region_ids' in inputs:
        region_ids = [int(r) for r in region_ids.split(',')]
    else:
        region_ids = np.unique(regions[~nodata_mask])
    
    ar_row, ar_col = np.indices(regions.shape, dtype=np.int32)
    if test_txt:
        test_sample = pd.read_csv(test_txt)
        if 'row' not in test_sample.columns:
            test_sample = pd.read_csv(test_txt, sep='\t')
    else:  
        train_sample = pd.read_csv(train_txt, sep='\t', index_col='obs_id')
        
        # Set any pixels used for training to -1 so they can be avoided for testing
        
        ar_row[train_sample.row, train_sample.col] = -1
        ar_col[train_sample.row, train_sample.col] = -1
        ar_row[nodata_mask] = -1
        ar_col[nodata_mask] = -1
        
        test_rows = random.sample(ar_row[ar_row > -1], n_samples)
        test_cols = random.sample(ar_col[ar_col > -1], n_samples)
        test_sample = pd.DataFrame({'row': test_rows, 'col': test_cols})
        test_sample['region'] = regions[test_rows, test_cols]
        test_basename = 'test_sample_%s.txt' % n_samples
        test_sample.to_csv(os.path.join(out_dir, test_basename))
        desc = 'Random test sample of %s not used in training samples %s' % (n_samples, train_txt)
        createMetadata(sys.argv, os.path.join(out_dir, test_basename), description=desc)
    
    ind_mask = ar_row != -1
    
    confusion_dir = os.path.join(out_dir, 'region_confusion_tbls')
    if not os.path.isdir(os.path.join(confusion_dir)):
        os.mkdir(confusion_dir)
        
    stats = []
    n_regions = len(region_ids)
    for i, r_id in enumerate(region_ids):
        print '\nCalculating stats for %s of %s regions...' % (i + 1, n_regions)
        t1 = time.time()
        region_mask = regions == r_id
        min_row = int(ar_row[region_mask & ind_mask].min())
        max_row = int(ar_row[region_mask & ind_mask].max())
        min_col = int(ar_col[region_mask & ind_mask].min())
        max_col = int(ar_col[region_mask & ind_mask].max())
        nrows = max_row - min_row
        ncols = max_col - min_col

        ar_t = ds_t.ReadAsArray(min_col, min_row, ncols, nrows)
        ar_p = ds_p.ReadAsArray(min_col, min_row, ncols, nrows)
        clipped_mask = region_mask[min_row:max_row, min_col:max_col] 
        clipped_mask = (ar_t == t_nodata) | (ar_p == p_nodata)
        del region_mask
        
        region_sample = test_sample[test_sample.region == r_id].copy()
        region_sample['global_row'] = region_sample.row
        region_sample['global_col'] = region_sample.col
        region_sample['row'] = region_sample.row - min_row
        region_sample['col'] = region_sample.col - min_col
        region_sample['reference'] = ar_t[region_sample.row, region_sample.col]
        region_sample['predicted'] = ar_p[region_sample.row, region_sample.col]
        
        
        df = evaluation.confusion_matrix_by_area(ar_p, ar_t, region_sample, p_nodata, t_nodata, mask=clipped_mask, bins=bins, match='best')
        this_txt = os.path.join(confusion_dir, 'confusion_%s.txt' % r_id)
        df.to_csv(this_txt, sep='\t')
        accuracy = df.ix['producer', 'user']
        kappa = df.ix['producer', 'kappa']
        rmse = sf.rmse(region_sample.reference, region_sample.predicted)
        
        stats.append({'region': r_id,
                      'accuracy': accuracy,
                      'kappa': kappa,
                      'rmse': rmse
                      })
        print 'Time for this region: %.1f minutes' % ((time.time() - t1)/60)
    df = pd.DataFrame(stats)
    out_txt = os.path.join(out_dir, 'region_stats.txt')
    df.to_csv(out_txt, sep='\t', index=False)
    test_basename = os.path.basename(test_txt)
    desc = 'Stats from modeling regions from %s. Samples drawn from %s' % (region_raster, os.path.join(out_dir, test_basename))
    createMetadata(sys.argv, out_txt, description=desc)
    
    if not test_txt:
        print '\nTest sample text file written to:', test_txt

    if region_shp:
        #df['FID'] = df.region
        #import pdb; pdb.set_trace()
        out_shp = os.path.join(out_dir, 'modeling_region_stats.shp')
        df_to_shp(df, region_shp, out_shp, copy_fields=False, df_id_field='region', shp_id_field='Zone_ID')
    print '\nFinished in %.1f minutes' % ((time.time() -t0)/60)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
        
        
        
        
