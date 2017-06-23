# -*- coding: utf-8 -*-
"""
Created on Fri May 26 23:52:24 2017

@author: shooper
"""
import gdal
import numpy as np
import seaborn as sns

from evaluation import get_samples, histogram_2d

def main():
    p_path = '/vol/v2/stem/imperv/models/imperv_20161012_0958/imperv_20161012_0958_vote.bsq'
    t_path = '/vol/v2/stem/imperv/truth_map/imperv2001_CAORWA.bsq'
    
    nodata_p = 255
    nodata_t = 255
    out_dir = '/vol/v2/stem/imperv/models/imperv_20161012_0958/evaluation_vote/'
    
    ds_p = gdal.Open(p_path)
    ar_p = ds_p.ReadAsArray()
    ds_p = None
    
    ds_t = gdal.Open(t_path)
    ar_t = ds_t.ReadAsArray()
    ds_t = None
    
    sample_txt = '/vol/v2/stem/imperv/samples/imperv_sample1454542_20161007_0843/imperv_sample1454542_20161007_0843_test.txt'
    #sample_txt = '/vol/v2/stem/canopy/samples/canopy_sample1454542_20161017_1919/canopy_sample1454542_20161017_1919_test.txt'
    df = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    p_samples, t_samples = get_samples(ar_p, ar_t, df, nodata_p, nodata_t, match=False)
    #t_mask = (t_samples > 0) & (p_samples > 0)
    #t_samples = t_samples[t_mask]
    #p_samples = p_samples[t_mask]
    out_png = os.path.join(out_dir, 'imperv_20161012_0958_2dhistogram_average_hex_gray.png')
    #out_png = os.path.join(out_dir, 'canopy_20161018_2254_2dhistogram_bestmatch_hex_gray.png')
    #import pdb; pdb.set_trace()
    
    sns.set_context(context='paper', font_scale=1.4)
    histogram_2d(t_samples, p_samples, out_png, bins=50, hexplot=True, vmax=4000)
    print out_png
    ar_p = None
    ar_t = None
    p_samples = None
    t_samples = None

main()