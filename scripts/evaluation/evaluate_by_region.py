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
    #nodata_mask = regions == region_nodata
    ds = None
    
    ds_t = gdal.Open(t_raster)
    ar_t = ds_t.ReadAsArray()
    
    ds_p = gdal.Open(p_raster)
    ar_p = ds_p.ReadAsArray()    
    nodata_mask = (ar_p == p_nodata) | (ar_t == t_nodata)
    
    d = {}
    if 'region_ids' in inputs:
        region_ids = [int(r) for r in region_ids.split(',')]
    else:
        region_ids = np.unique(regions[~nodata_mask])
    
    if test_txt:
        test_sample = pd.read_csv(test_txt)
        if 'row' not in test_sample.columns:
            test_sample = pd.read_csv(test_txt, sep='\t')
    else:  
        train_sample = pd.read_csv(train_txt, sep='\t', index_col='obs_id')
        
        # Set any pixels used for training to -1 so they can be avoided for testing
        n_rows, n_cols = regions.shape
        n_pixels = regions.size
        pixel_ids = np.arange(n_pixels, dtype=np.uint32).reshape(n_rows, n_cols)
        pixel_ids[train_sample.row, train_sample.col] = n_pixels #will always be 1 more than last col
        pixel_ids[nodata_mask] = n_pixels
        #ar_col[train_sample.row, train_sample.col] = -1

        #import pdb; pdb.set_trace()
        
        
        test_ids = np.array(random.sample(pixel_ids[pixel_ids != n_pixels], n_samples), dtype=np.uint32)
        test_rows = test_ids / n_cols
        test_cols = test_ids % n_cols
        #test_cols = random.sample(ar_col[ar_col != -1], n_samples)
        test_sample = pd.DataFrame({'row': test_rows, 'col': test_cols})
        test_sample['region'] = regions[test_rows, test_cols]
        
    
    ind_mask = (pixel_ids == n_pixels).reshape(n_rows, n_cols)
    del pixel_ids
    ar_row, ar_col = np.indices(regions.shape, dtype=np.int32)
    
    confusion_dir = os.path.join(out_dir, 'region_confusion_tbls')
    if not os.path.isdir(os.path.join(confusion_dir)):
        os.mkdir(confusion_dir)
    
    test_sample['reference'] = -1
    test_sample['predicted'] = -1
    stats = []
    n_regions = len(region_ids)
    for i, r_id in enumerate(region_ids):
        print '\nCalculating stats for region %s (%s of %s)...' % (r_id, i + 1, n_regions)
        t1 = time.time()
        ''' figure out what's up with nodata values '''
        region_mask = regions == r_id
        min_row = int(ar_row[region_mask & ~ind_mask].min())
        max_row = int(ar_row[region_mask & ~ind_mask].max())
        min_col = int(ar_col[region_mask & ~ind_mask].min())
        max_col = int(ar_col[region_mask & ~ind_mask].max())
        nrows = max_row - min_row
        ncols = max_col - min_col

        ar_t_region = ds_t.ReadAsArray(min_col, min_row, ncols, nrows)
        ar_p_region = ds_p.ReadAsArray(min_col, min_row, ncols, nrows)
        #clipped_mask = region_mask[min_row:max_row, min_col:max_col] 
        clipped_mask = (ar_t_region == t_nodata) | (ar_p_region == p_nodata)
        del region_mask
        
        region_sample = test_sample[test_sample.region == r_id].copy()
        region_sample['global_row'] = region_sample.row
        region_sample['global_col'] = region_sample.col
        region_sample['row'] = region_sample.row - min_row
        region_sample['col'] = region_sample.col - min_col
        region_sample['reference'] = ar_t_region[region_sample.row, region_sample.col]
        region_sample['predicted'] = ar_p_region[region_sample.row, region_sample.col]
    
        region_sample = region_sample
        test_sample.ix[region_sample.index, 'reference'] = region_sample['reference']
        test_sample.ix[region_sample.index, 'predicted'] = region_sample['predicted']
        #import pdb; pdb.set_trace()
        #try:
        df = evaluation.confusion_matrix_by_area(ar_p_region, ar_t_region, region_sample, p_nodata, t_nodata, mask=clipped_mask, bins=bins, match='best')
        this_txt = os.path.join(confusion_dir, 'confusion_%s.txt' % r_id)
        df.to_csv(this_txt, sep='\t')
        accuracy = df.ix['producer', 'user']
        kappa = df.ix['producer', 'kappa']
        sample_mask = (region_sample.reference==t_nodata) | (region_sample.predicted==p_nodata)
        rmse = sf.rmse(region_sample.reference[~sample_mask], region_sample.predicted[~sample_mask])
        print len(sample_mask[sample_mask])
        
        stats.append({'region': r_id,
                      'accuracy': accuracy,
                      'kappa': kappa,
                      'rmse': rmse
                      })
        print 'Time for this region: %.1f minutes' % ((time.time() - t1)/60)

    #if not test_txt:
    test_basename = 'test_sample_%s.txt' % n_samples
    test_txt = os.path.join(out_dir, test_basename)
    test_sample.to_csv(test_txt)
    desc = 'Random test sample of %s not used in training samples %s' % (n_samples, train_txt)
    createMetadata(sys.argv, os.path.join(out_dir, test_basename), description=desc)
    print '\nTest sample text file written to:', test_txt    
    
    df = pd.DataFrame(stats)
    out_txt = os.path.join(out_dir, 'region_stats_%s.txt' % n_samples)
    df.to_csv(out_txt, sep='\t', index=False)
    test_basename = os.path.basename(test_txt)
    desc = 'Stats from modeling regions from %s. Samples drawn from %s' % (region_raster, os.path.join(out_dir, test_basename))
    createMetadata(sys.argv, out_txt, description=desc)

    if region_shp:
        out_shp = os.path.join(out_dir, 'modeling_region_stats.shp')
        df_to_shp(df, region_shp, out_shp, copy_fields=False, df_id_field='region', shp_id_field='Zone_ID')
    print '\nFinished in %.1f minutes' % ((time.time() -t0)/60)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
        
        
        
        
