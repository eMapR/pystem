# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 16:16:08 2016

@author: shooper
"""

import os
import gdal
import fnmatch
import glob
import pandas as pd
import numpy as np

import stem


def main():
    srch_dir = '/vol/v2/stem/canopy/models'
    stamps = fnmatch.filter(os.listdir(srch_dir), 'canopy*')
    
    info = []
    for stamp in stamps:
        print stamp
        this_dir = os.path.join(srch_dir, stamp)
        this_srch_str = os.path.join(this_dir, 'train_stem*_params.txt')
        matched = glob.glob(this_srch_str)
        if len(matched) == 0:
            print 'No param file for ', stamp
            info.append([stamp,
                         '',
                         0,
                         0,
                         0,
                         0,
                         False,
                         0,
                         0,
                         False,
                         0,
                         0,
                         '',
                         0,
                         ''
                         ])
            continue
        
        this_param_text = matched[0]
        if 'regressor' in this_param_text: 
            model_type = 'Regressor'
        else: 
            model_type = 'Classifier'
        inputs, df_var = stem.read_params(this_param_text)
        for var in inputs:
            exec ("{0} = str({1})").format(var, inputs[var])
        
        vote_mask = False
        vote_accuracy = None
        vote_kappa = None
        vote_dir = os.path.join(this_dir, 'evaluation_vote')
        if os.path.exists(vote_dir):
            for root, dis, files in os.walk(vote_dir):
                for f in files:
                    if f.endswith('txt'):
                        vote_txt = os.path.join(root, f)
                        df_vote = pd.read_csv(vote_txt, sep='\t', index_col='bin')
                        try:
                            vote_accuracy = int(df_vote.ix['producer','user'])
                            vote_kappa = round(df_vote.ix['kappa','kappa'], 2)
                        except:
                            vote_accuracy = int(df_vote.ix['user','producer'])
                            vote_kappa = round(df_vote.ix['user','kappa'], 2)                            
                        if 'mask' in vote_txt: vote_mask=True
        
        mean_mask = False
        mean_accuracy = None
        mean_kappa = None
        mean_dir = os.path.join(this_dir, 'evaluation_mean')
        if os.path.exists(mean_dir):
            for root, dis, files in os.walk(mean_dir):
                for f in files:
                    if f.endswith('txt'):
                        mean_txt = os.path.join(root, f)
                        df_mean = pd.read_csv(mean_txt, sep='\t', index_col='bin')
                        try:
                            mean_accuracy = int(df_mean.ix['producer','user'])
                            mean_kappa = round(df_mean.ix['kappa','kappa'], 2)
                        except:
                            mean_accuracy = int(df_mean.ix['user','producer'])
                            mean_kappa = round(df_mean.ix['user','kappa'], 2)                           
                        if 'mask' in mean_txt: mean_mask=True
        
        dt_dir = os.path.join(this_dir, 'decisiontree_models')
        try:
            n_sets = len(os.listdir(dt_dir)) - 1
        except:
            n_sets = None
        
        n_samples = int(sample_txt.split('_')[1].replace('sample',''))
        
        if not 'max_features' in inputs: max_features = None
        
        
        avg_count = None
        cnt_path = os.path.join(this_dir, '%s_count.bsq' % stamp)
        if os.path.exists(cnt_path):
            ds = gdal.Open(cnt_path)
            ar = ds.ReadAsArray()
            cnt_nodata = ds.GetRasterBand(1).GetNoDataValue()
            ds = None
            if len(ar[ar == cnt_nodata]) == 0:
                cnt_min = ar.min()
                cnt_max = ar.max()
                if cnt_min <= 0: 
                    cnt_nodata = cnt_min
                else: 
                    cnt_nodata = cnt_max
            avg_count = int(round(np.mean(ar[ar != cnt_nodata]),0))   
        
        avg_oob = None
        oob_path = os.path.join(this_dir, '%s_oob.bsq' % stamp)
        if os.path.exists(oob_path):
            ds = gdal.Open(oob_path)
            ar = ds.ReadAsArray()
            ds = None
            oob_min = ar.min()
            oob_max = ar.max()
            if oob_min <= 0: 
                oob_nodata = oob_min
            else: 
                oob_nodata = oob_max
            avg_oob = round(np.mean(ar[ar != oob_nodata]),1)
            

        
        info.append([stamp,
                     model_type,
                     avg_oob,
                     avg_count,
                     vote_accuracy,
                     vote_kappa,
                     vote_mask,
                     mean_accuracy,
                     mean_kappa,
                     mean_mask,
                     n_sets,
                     n_samples,
                     '[%s]' % support_size,
                     sets_per_cell,
                     max_features
                     ])
                     
    df = pd.DataFrame(info, columns=['stamp',
                                     'model_type',
                                     'avg_oob',
                                     'avg_count',
                                     'vote_accuracy', 
                                     'vote_kappa', 
                                     'vote_mask', 
                                     'mean_accuracy', 
                                     'mean_kappa',
                                     'mean_mask', 
                                     'n_sets', 
                                     'n_samples',
                                     'support_size',
                                     'sets_per_cell',
                                     'max_features'])
    
    out_txt = os.path.join(srch_dir, 'model_info.txt')
    df.to_csv(out_txt, sep='\t', index=False)
    
    print 'Text written to ', out_txt

main()
        
        
        
    