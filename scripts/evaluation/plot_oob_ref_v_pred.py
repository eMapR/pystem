# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:42:50 2017

@author: shooper
"""

import os
import cPickle as pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random

import stem

from lthacks import attributes_to_df
from lthacks import stats_functions as sf




def main(set_txt, train_params, plot_dims, out_png):
    
    sns.set_context(context='paper', font_scale=.3)
    sns.set_style('white', rc={'axes.linewidth': .5})
    
    inputs, df_var = stem.read_params(train_params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])     
    predict_cols = sorted(df_var.index)
    df = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    df = df[df.oob_rate < 5]

    fig, axes = plt.subplots(*plot_dims)
    n_plots = axes.size
    set_ids = random.sample(df.index, n_plots)
    
    sample = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    
    for si, ax in zip(set_ids, axes.ravel()):
        support_set = df.ix[si]
        
        oob_ind_txt = support_set.dt_file.replace('.pkl', '_oob_inds.txt')
        with open(oob_ind_txt) as txt: 
            oob_inds = [int(l) for l in txt]
        with open(support_set.dt_file, 'rb') as f:
            dt_model = pickle.load(f)
        
        oob_sample = sample.ix[oob_inds]
        oob_predictions = dt_model.predict(oob_sample[predict_cols])
        oob_reference = oob_sample[target_col]
        
        rmse = sf.rmse(oob_reference, oob_predictions)
        ac, ac_s, ac_u, ssd, spod = sf.agree_coef(oob_reference, oob_predictions)
        ax.plot(oob_reference, oob_predictions, 'o', alpha=0.05, markeredgecolor='none', markersize=2.5)
        #ax.xticks([0, 50, 100])
        #ax.yticks([0, 50, 100])
        title = 'Set ID: %s, RMSE: %.1f, ac: %.3f' % (si, rmse, ac)
        ax.set_title(title)
        sns.despine()
    
    fig.subplots_adjust(hspace=0.1)
    plt.savefig(out_png, dpi=300)
        
        
        
set_txt = '/vol/v2/stem/conus/models/imperv_20171025_1149/decisiontree_models/imperv_20171025_1149_support_sets.txt'
train_params =  '/vol/v2/stem/conus/models/imperv_20171025_1149/train_stem_params.txt'
out_png = os.path.join(os.path.dirname(train_params), 'low_oob_plot.png')
plot_dims = 5, 5
main(set_txt, train_params, plot_dims, out_png)    