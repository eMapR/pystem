# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:33:36 2017

@author: shooper
"""

import os
import sys
import gdal
import glob
import numpy as np
import pandas as pd
from confusion_matrix import read_params
import stem
from evaluation import calc_rmse, get_samples


def main(search_dir, models, t_path, inventory_txt):
    
    df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
    columns = df_inv.columns
    if 'vote_rmse' not in columns:
        df_inv['vote_rmse'] = None
    if 'mean_rmse' not in columns:
        df_inv['mean_rmse'] = None
    df_inv = df_inv.ix[models]
    
    ds = gdal.Open(t_path)
    ar_t = ds.ReadAsArray()
    ds = None

    for model in models:
        print '\nCalculating RMSE for ', model
        model_dir = os.path.join(search_dir, model)
        if not os.path.exists(model_dir):
            print 'Model dir does not exist: %s. Skipping...\n' % model_dir
            continue
        
        confusion_params = os.path.join(model_dir, 'confusion_matrix_params.txt')
        if not os.path.exists(confusion_params):
            print 'Could not find confusion params: ', confusion_params
            predict_params = os.path.join(model_dir, 'predict_stem_params.txt')
            inputs, _ = stem.read_params(predict_params)
            p_nodata = int(inputs['nodata'].replace('"',''))
            this_srch_str = os.path.join(model_dir, 'train_stem*_params.txt')
            train_params = glob.glob(this_srch_str)
            if len(train_params) == 0: 
                print 'Can not find test data for ', model, '\n'
                continue
            train_params = train_params[0]
            inputs, _ = stem.read_params(train_params)
            test_txt = inputs['sample_txt'].replace('predictors', 'test').replace('"','')
        else:
            inputs = read_params(confusion_params)
            for k, v in inputs.iteritems():
                inputs[k] = v.replace('"','')
            test_txt = inputs['sample_txt']
            p_nodata = int(inputs['p_nodata'])
        df = pd.read_csv(test_txt, sep='\t', index_col='obs_id')
        for agg_method in ['vote', 'mean']:
            p_path = os.path.join(model_dir, '%s_%s.bsq' % (model, agg_method))
            ds = gdal.Open(p_path)
            ar_p = ds.ReadAsArray()
            t_samples, p_samples = get_samples(ar_p, ar_t, p_nodata, 255, samples=df, match=True)
            #import pdb; pdb.set_trace()
            
            df_inv.ix[model, '%s_rmse' % agg_method] = np.round(calc_rmse(t_samples, p_samples),1)

    out_txt = inventory_txt.replace('.txt', '_final.txt')
    df_inv.to_csv(out_txt, sep='\t')

'''search_dir = '/vol/v2/stem/canopy/models'
models = [
'canopy_20161009_1058', 
'canopy_20161009_2111', 
'canopy_20161017_1358',
'canopy_20161013_1435',
'canopy_20161017_2013',
'canopy_20161018_2254',
'canopy_20161111_1413']
t_path = '/vol/v2/stem/canopy/truth_map/canopy2001_CAORWA.bsq'
inventory_txt = os.path.join(search_dir, 'model_info.txt')
main(search_dir, models, t_path, inventory_txt)#'''


search_dir = '/vol/v2/stem/imperv/models'
models = [
'imperv_20161019_0823',
'imperv_20161018_2307',
'imperv_20160929_1617',
'imperv_20161002_1906',
'imperv_20161007_2240',
'imperv_20161012_0958',
'imperv_20161012_2350'
]

#models = ['imperv_20161007_2240']
t_path = '/vol/v2/stem/imperv/truth_map/imperv2001_CAORWA.bsq'
inventory_txt = os.path.join(search_dir, 'model_info.txt')
main(search_dir, models, t_path, inventory_txt)#'''