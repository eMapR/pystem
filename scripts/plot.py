# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 21:09:53 2016

@author: SamHooper
"""
import os
import gdal
import glob
import fnmatch
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


def plot_per_row(df, out_png, index_col=None, xlabels=None, label_col=None, xlabel='Evaluation Metric', ylabel='Score', title=None, ylim=None, bar_width=.25):
    print 'Plotting stats...'
    #plt.style.use('ggplot')
    sns.set_style('white')
    sns.set_context(context='talk', font_scale=1.4, rc={'patch.linewidth': 0, 'face.color':79/255.0})
    
    if index_col:
        df.set_index(index_col, inplace=True)
    columns = [c for c in df.columns if c not in ['color', 'hatch']]
    n_cols = len(columns)
    n_rows = len(df)
    class_width = n_rows * bar_width
    x_loc = np.arange(n_cols) * (class_width + bar_width) + class_width * .75
    if n_rows % 2 == 0:
        df['placement'] = range(-1 * n_rows/2, n_rows/2)
    else:
        df['placement'] = range(-1 * (n_rows/2), n_rows/2 + 1)

    for i, r in df.iterrows():
        if 'label' in df.columns:
            label = r.label
        elif label_col:
            label = r[label_col]
        else:
            label = i
        hatch = None
        edgecolor=None
        if 'hatch' in df.columns: 
            hatch = r.hatch
            edgecolor='k'
        bar = plt.bar(x_loc + (bar_width * r.placement), r[columns],  bar_width, color=r.color, ecolor='none', label=label, hatch=hatch)


    #ax = plt.gca()
    #ax.set_axis_bgcolor('.31')
    #ax.axes.get_yaxis().set_visible(False)
    #ax.axes._frameon = False
    #ax.axes.get_xaxis().set_visible(False)

    if xlabels:
        plt.xticks(x_loc + bar_width / 2. - bar_width * 2, xlabels, size='small', rotation=45, multialignment='right')
        plt.subplots_adjust(bottom=0.2)
        plt.tick_params(axis='x', width=0)#, pad=-2)
    else:
        plt.tick_params(axis='x', which='minor', color='none')
        #plt.xticks(x_loc + width / 2., size='small')

    plt.tick_params(labelsize=10, right=False)
    plt.tick_params(axis='y', which='minor', color='none')
    #plt.ylabel(ylabel)
    #plt.xlabel(xlabel, labelpad=-3)

    if title:
        plt.title(title, color='.3', weight='normal', fontsize=12)
    #plt.xlim(0, max(x_loc) + class_width/2 + bar_width * 8)
    plt.xlim(0, max(x_loc) + class_width/2)# + bar_width * 4)

    if ylim:
        plt.ylim(ylim[0], ylim[1])

    #plt.plot((0, max(x_loc) + width * 2), (0, 0), color='k')
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    legend = plt.legend(loc='lower left', bbox_to_anchor=(1, -0.05), frameon=False)
    for text in legend.get_texts():
        plt.setp(text, color='.3')#'''
    #plt.show()
    sns.despine()
    sns.axes_style({'axes.linewidth': 0.4})
    
    fig = plt.gcf()
    #fig.set_size_inches(4.1, 3.925)
    plt.savefig(out_png, dpi=300)
    plt.clf()
    print 'Plot saved to ', out_png


def plot_by_column(df, out_png, color_dict, new_index, labels=None, xlabel='Evaluation Metric',
                   ylabel='Score', title=None, ylim=None, bar_width=.25):
    print 'Plotting stats...'
    plt.style.use('ggplot')

    df = df.set_index(new_index)
    n_cols = len(df.columns)
    class_width = n_cols * bar_width
    x_loc = np.arange(n_cols) + class_width / 2.0
    for i in df.columns():

        plt.bar()
    bar_05 = plt.bar(x_loc - width, df['5'],  width, color=color_dict['5'],  edgecolor='none')
    bar_30 = plt.bar(x_loc,         df['30'], width, color=color_dict['30'], edgecolor='none', label='30 m')
    bar_60 = plt.bar(x_loc + width, df['60'], width, color=color_dict['60'], edgecolor='none', label='60 m')

    #if labels != None:
    plt.xticks(x_loc + width / 2. - width, labels, size='small', rotation=45,multialignment='right')
    plt.subplots_adjust(bottom=0.2)
    #else:
        #plt.xticks(x_loc + width / 2., size='small')

    plt.tick_params(labelsize=8, right=False)
    plt.tick_params(axis='x', width=0, pad=-7)
    plt.tick_params(axis='y', which='minor', color='none', width=0)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    if title:
        plt.title(title, color='.3', weight='normal', fontsize=12)
    plt.xlim(0, max(x_loc) + width * 2)
    if ylim:
        plt.ylim(ylim[0], ylim[1])

    #plt.plot((0, max(x_loc) + width * 2), (0, 0), color='k')

    legend = plt.legend(loc='upper right', frameon=False)
    for text in legend.get_texts():
        plt.setp(text, color='.3', labelsize=8)
    # plt.legend(bar_s[0], bar_u[0], ['Systematic', 'Unsystematic'])

    plt.savefig(out_png)
    plt.clf()
    print 'Plot saved to ', out_png


#def plot_auc_by_lc(txt1, txt2, txt5, txt10, out_png, color_dict=None):
def plot_stats_by_lc(search_dir, search_str, out_png, color_dict=None):

    sns.set_style('white')
    sns.set_context('paper', font_scale=2.5)

    files = fnmatch.filter(os.listdir(search_dir), search_str)
    for f in files:
        if 'res1yr' in f:
            txt1 = os.path.join(search_dir, f)
        if 'res2yr' in f:
            txt2 = os.path.join(search_dir, f)
        if 'res5yr' in f:
            txt5 = os.path.join(search_dir, f)
        if 'res10yr' in f:
            txt10 = os.path.join(search_dir, f)
    
    df1 = pd.read_csv(txt1, sep='\t', index_col='lc_class')
    df2 = pd.read_csv(txt2, sep='\t', index_col='lc_class')
    df5 = pd.read_csv(txt5, sep='\t', index_col='lc_class')
    df10 = pd.read_csv(txt10, sep='\t', index_col='lc_class')
    
    pnl = pd.Panel({1:df1, 2:df2, 5:df5, 10:df10})
    df_rmse = pnl[:, :, 'rmse']
    df_auc = pnl[:, :, 'auc']
    
    label_dict = {23: 'urban',
                  11: 'water',
                  31: 'bare ground',
                  41: 'deciduous',
                  42: 'coniferous',
                  52: 'shrub'}
    
    for lc in df1.index:
        if color_dict:
            plt.plot(df_auc.columns, 
                     df_auc.ix[lc], 
                     color=color_dict[lc], 
                     label=label_dict[lc])
        else:
            plt.plot(df_auc.columns, df_auc.ix[lc], label=label_dict[lc])
        plt.xticks(df_auc.columns)
        #plt.yticks([0, 0.25, 0.5, 0.75, 1.0])
        plt.xlim(0, 11)
        plt.ylim(0, 1)
        #plt.legend()
    sns.despine()
    plt.savefig(out_png.replace('.png', '_auc.png'))
    plt.clf()
    
    for lc in df1.index:
        if color_dict:
            plt.plot(df_rmse.columns, 
                     1 - df_rmse.ix[lc]/100.0, 
                     color=color_dict[lc], 
                     label=label_dict[lc])     
        else:
            plt.plot(df_rmse.columns, 1 - df_rmse.ix[lc]/100.0, label=label_dict[lc])
        plt.xticks(df_rmse.columns)
        plt.xlim(0, 11)
        plt.ylim(0, 1)
        #plt.legend()
    sns.despine()
    sns.axes_style({'axes.linewidth': .7})
    plt.savefig(out_png.replace('.png', '_rmse.png'))
    print 'Plots saved to: ', os.path.dirname(out_png)


def obs_per_lc(txt, lc_path, target_col, out_png, xlabels=None):
    
    ds = gdal.Open(lc_path)
    ar = ds.ReadAsArray()
    ds = None

    df = pd.read_csv(txt, sep='\t', index_col='obs_id')
    df = df[df.YEAR == 2011]
    df['lc'] = ar[df.row, df.col]
    lc_vals = np.array([11,12,23,31,41,42,52])
    bins = sorted(lc_vals) + [lc_vals.max() + 1]
    c1, b = np.histogram(df.ix[df[target_col] == 1, 'lc'], bins=bins)
    c0, b = np.histogram(df.ix[df[target_col] == 0, 'lc'], bins=bins)

    for lc in bins[:-1]:
        print 'land cover: ', lc
        print 'positive: ', len(df[(df.lc == lc) & (df[target_col] == 1)])
        print 'negative: ', len(df[(df.lc == lc) & (df[target_col] == 0)])
        print ''
    plt.clf()
    loc = np.arange(1, len(c0) + 1)
    plt.bar(loc, c0, width=.35, color='k', label='neg')
    plt.bar(loc, c1, width=.35, bottom=c0, color='.5', label='pos')
    if not xlabels:
        xlabels = bins
    plt.xticks(loc + .35/2, xlabels, size='small', rotation=45, multialignment='right')
    n = len(df)
    title = '%s,  n = %s' % (os.path.basename(os.path.dirname(txt)).split('_res')[0], n)
    plt.title(title)
    plt.legend()
    plt.savefig(out_png)


def plot_importance(model_txt, species_list, search_dir, constant_vars, out_png):
    
    sns.set_style('white')
    sns.set_context('talk')
    
    df = pd.read_csv(model_txt, index_col='stamp')
    #models = [model for model in df.index if species in model]
    #df = df.ix[models]
    df['constant'] = 0
    df['yearly'] = 0
    df['obs_year'] = 0
    
    for model in df.index:
        this_dir = os.path.join(search_dir, model)
        var_txt = glob.glob(os.path.join(this_dir, 'var_text_ebird*.txt'))[0]
        df_var = pd.read_csv(var_txt, sep='\t', index_col='var_name')
        constant =  df_var.ix[constant_vars, 'importance'].mean()
        yearly = df_var.ix[~df_var.index.isin(constant_vars), 'importance'].mean()
        df.ix[model, 'constant'] = constant
        df.ix[model, 'yearly'] = yearly
        df.ix[model, 'obs_year'] = df_var.ix['YEAR', 'importance']
    
    for species in species_list:
        models = df.ix[[model for model in df.index if species in model]]
        plt.plot(models.temporal_res, 
                 models.constant, 
                 '-', 
                 label='constant',
                 alpha=0.4,
                 color=[0.769, 0.306, 0.322])
                 
        plt.plot(models.temporal_res,
                 models.yearly,
                 '-',
                 label='yearly',
                 alpha=.4,
                 color=[0.298, 0.451, 0.69])
        plt.plot(models.temporal_res,
                 models.obs_year,
                 alpha=.4,
                 color='k')
        
        #'''
    resolutions = df.temporal_res.unique()
    """summary = []
    for res in resolutions:
        df_res = df[df.temporal_res == res]
        summary.append({'constant_mean': df_res.constant.mean(),
                        'constant_max': df_res.constant.max() - df_res.constant.mean(),
                        'constant_min': df_res.constant.mean() - df_res.constant.min(),
                        'constant_stdv': df_res.constant.std(),
                        'yearly_mean': df_res.yearly.mean(),
                        'yearly_max': df_res.yearly.max() - df_res.yearly.mean(),
                        'yearly_min': df_res.yearly.mean() - df_res.yearly.min(),
                        'yearly_stdv': df_res.yearly.std(),
                        'obs_yr_mean': df_res.obs_year.mean(),
                        'obs_yr_max': df_res.obs_year.max() - df_res.obs_year.mean(),
                        'obs_yr_min': df_res.obs_year.mean() - df_res.obs_year.min(),
                        'obs_yr_stdv': df_res.obs_year.std()})
        '''constant_mean = df_res.constant.mean()
        constant_max = df_res.constant.max() - constant_mean
        constant_min = constant_mean - df_res.constant.min()
        constant_stdv = df_res.constant.std()
        yearly_mean = df_res.yearly.mean() 
        yearly_max = df_res.yearly.max() - yearly_mean
        yearly_min = yearly_mean - df_res.yearly.min()
        constant_stdv = df_res.consta
        obs_yr_mean = df_res.obs_year.mean()
        obs_yr_max = df_res.obs_year.max() - obs_yr_mean
        obs_yr_min = obs_yr_mean - df_res.obs_year.min()
        summary.append({'constant_mean': constant_mean,
                        'constant_max': constant_max,
                        'constant_min': constant_min,
                        'yearly_mean': yearly_mean,
                        'yearly_max': yearly_max,
                        'yearly_min': yearly_min,
                        'obs_yr_mean': obs_yr_mean,
                        'obs_yr_max': obs_yr_max,
                        'obs_yr_min': obs_yr_min})'''
    df_summary = pd.DataFrame(summary)

    '''plt.plot(resolutions,
                 df_summary.constant_mean,
                 #yerr=[df_summary.constant_min, df_summary.constant_max],
                 yerr=df_summary.constant_stdv,
                 alpha=0.3,
                 fmt='o',
                 capthick=1,
                 ecolor=[0.769, 0.306, 0.322],
                 color=[0.769, 0.306, 0.322])'''
    plt.plot(resolutions,
             df_summary.constant_mean,
             label='constant',
             color=[0.769, 0.306, 0.322])
    
    '''plt.errorbar(resolutions,
                 df_summary.yearly_mean,
                 #yerr=[df_summary.yearly_min, df_summary.yearly_max],
                 yerr=df_summary.yearly_stdv,
                 alpha=0.3,
                 fmt='o',
                 capthick=1,
                 ecolor=[0.298, 0.451, 0.69],
                 color=[0.298, 0.451, 0.69])'''
    plt.plot(resolutions,
             df_summary.yearly_mean,
             label='yearly',
             color=[0.298, 0.451, 0.69])

    '''plt.errorbar(resolutions,
                 df_summary.obs_yr_mean,
                 #yerr=[df_summary.obs_yr_min, df_summary.obs_yr_max],
                 yerr=df_summary.obs_yr_stdv,
                 alpha=0.3,
                 fmt='o',
                 capthick=1,
                 ecolor=[0.3, 0.3, 0.3],
                 color=[0.3, 0.3, 0.3])'''
    plt.plot(resolutions,
             df_summary.obs_yr_mean,
             label='yearly',
             color=[0.3, 0.3, 0.3])"""
    
    #plt.legend()
    plt.xticks(resolutions, resolutions)
    plt.xlim(0, 11)
    #plt.ylim(0, .21)
    sns.despine()
    sns.axes_style({'axes.linewidth': .7})
    plt.savefig(out_png)
    
        

# Plot canopy and imperv stem/bdt stuff
txt = '/home/server/student/homes/shooper/stem_ch1_figures/stem_bdt.txt'
df = pd.read_csv(txt, sep=',')
colors = [[0.5, 0.5, 0.5],
          [0.9, 0.9, 0.9],
          [0.3, 0.3, 0.3],
          [0.7, 0.7, 0.7]
          ]
          
#purple/green          
'''colors = ['#268b85',
          '#6ec7c3',
          '#a5518c',
          '#c997ba']#'''
          
#burnt sienna/green          
colors = ['#0b7478',
          '#16b5b8',
          '#d03821',
          '#fc7a6a'
          ]#'''

#muted burnt sienna/green
colors = ['#1d7477',
          '#2badb0',
          '#c64936',
          '#f2877a'
          ]
df['color'] = colors
#df['hatch'] = [None, None, '/', '/']
out_png = '/home/server/student/homes/shooper/stem_ch1_figures/stem_vs_bdt_plot_color.png'
xlabels = [s.strip() for s in 'Average OOB, Accuracy, Kappa, 1 - RMSE'.split(', ')]
plot_per_row(df, out_png, index_col='model', xlabels=xlabels, title=None, ylim=[.6, 1])


# eBird stuff    
'''for species in ['c_guttatus', 'h_rustica', 'p_melanocephalus', 's_neglecta', 'z_leucophrys']:
    txt = '/vol/v2/stem/ebird/results/%s_metrics.csv' % species
    #stamp = os.path.basename(txt).split('_metrics')[0]
    df = pd.read_csv(txt, index_col='stamp')
    colors = [[179,0,0],
              [227,74,51],
              [252,141,89],
              [253,204,138]
              ]
    df['color'] = [[c/255.0 for c in color] for color in colors]
    #df['color'] = [[]]
    out_png = '/vol/v2/stem/ebird/results/%s_temporal_res_plot.png' % species
    xlabels = 'OOB\nscore,AUC,1 - RMSE'.split(',')
    plot_per_row(df, out_png, index_col='temporal_res', xlabels=xlabels, title=None, ylim=[.4, 1])#'''

"""txt1 = '/vol/v2/stem/ebird/models/z_leucophrys_res1yr_20161203_2153/z_leucophrys_res1yr_20161203_2153_eval_2011_land_cover.txt'
txt2 = '/vol/v2/stem/ebird/models/z_leucophrys_res2yr_20161203_2318/z_leucophrys_res2yr_20161203_2318_eval_2011_land_cover.txt'
txt5 = '/vol/v2/stem/ebird/models/z_leucophrys_res5yr_20161203_2257/z_leucophrys_res5yr_20161203_2257_eval_2011_land_cover.txt'
txt10 = '/vol/v2/stem/ebird/models/z_leucophrys_res10yr_20161203_2302/z_leucophrys_res10yr_20161203_2302_eval_2011_land_cover.txt'
out_png = '/vol/v2/stem/ebird/results/z_leucophrys_land_cover_plot.png'#"""

"""txt1 = '/vol/v2/stem/ebird/models/h_rustica_res1yr_20161203_2339/h_rustica_res1yr_20161203_2339_eval_2011_land_cover.txt'
txt2 = '/vol/v2/stem/ebird/models/h_rustica_res2yr_20161204_0026/h_rustica_res2yr_20161204_0026_eval_2011_land_cover.txt'
txt5 = '/vol/v2/stem/ebird/models/h_rustica_res5yr_20161204_0014/h_rustica_res5yr_20161204_0014_eval_2011_land_cover.txt'
txt10 = '/vol/v2/stem/ebird/models/h_rustica_res10yr_20161204_0022/h_rustica_res10yr_20161204_0022_eval_2011_land_cover.txt'#"""

'''txt1= '/vol/v2/stem/ebird/models/p_melanocephalus_res1yr_20161201_1315/p_melanocephalus_res1yr_20161201_1315_eval_2011_land_cover.txt'
txt2 = '/vol/v2/stem/ebird/models/p_melanocephalus_res2yr_20161205_0839/p_melanocephalus_res2yr_20161205_0839_eval_2011_land_cover.txt'
txt5 = '/vol/v2/stem/ebird/models/p_melanocephalus_res5yr_20161203_0855/p_melanocephalus_res5yr_20161203_0855_eval_2011_land_cover.txt'
txt10 = '/vol/v2/stem/ebird/models/p_melanocephalus_res10yr_20161205_1103/p_melanocephalus_res10yr_20161205_1103_eval_2011_land_cover.txt' #'''

'''txt1= '/vol/v2/stem/ebird/models/c_guttatus_res1yr_20161204_0038/c_guttatus_res1yr_20161204_0038_eval_2011_land_cover.txt'
txt2 = '/vol/v2/stem/ebird/models/c_guttatus_res2yr_20161204_0049/c_guttatus_res2yr_20161204_0049_eval_2011_land_cover.txt'
txt5 = '/vol/v2/stem/ebird/models/c_guttatus_res5yr_20161204_0046/c_guttatus_res5yr_20161204_0046_eval_2011_land_cover.txt'
txt10 = '/vol/v2/stem/ebird/models/c_guttatus_res10yr_20161204_0047/c_guttatus_res10yr_20161204_0047_eval_2011_land_cover.txt'
'''

colors = [[23, [135,126,126]],
          [31, [199, 173, 135]],
          [41, [48,134,39]],
          [42, [9, 83, 30]],
          [52, [129,168,95]]
          ]

colors = [[23, [209, 94, 94]],
          [11, [50, 50, 200]],
          [31, [191, 138, 38]],
          [41, [75,138,33]],
          [42, [0, 115, 105]],
          [52, [168,203,69]]
          ]


'''for search_str in ['c_guttatus*', 'h_rustica*', 'p_melanocephalus*', 's_neglecta*', 'z_leucophrys*']:
    search_dir = '/vol/v2/stem/ebird/results/performance_by_lc'
    #search_str = 'z_leucophrys*'
    out_png = '/vol/v2/stem/ebird/results/performance_by_lc/%s_land_cover_plot.png' % search_str.replace('*', '')
    color_dict = {k: [c/255.0 for c in rgb] for k, rgb in colors}
    plot_stats_by_lc(search_dir, search_str, out_png, color_dict)
    #plot_auc_by_lc(txt1, txt2, txt5, txt10, out_png, color_dict)#'''



'''#txt = '/vol/v2/stem/ebird/samples/p_melanocephalus_res1yr_20161201_1231/erd_v5_Pheucticus_melanocephalus_predictors_test.txt'
#txt = '/vol/v2/stem/ebird/samples/h_rustica_res1yr_20161203_2324/erd_v5_Hirundo_rustica_predictors_test.txt'
#txt = '/vol/v2/stem/ebird/samples/z_leucophrys_res1yr_20161203_2147/erd_v5_Zonotrichia_leucophrys_predictors_test.txt'
#txt = '/vol/v2/stem/ebird/samples/c_guttatus_res1yr_20161204_0030/erd_v5_Catharus_guttatus_predictors_test.txt'
txt = '/vol/v2/stem/ebird/samples/s_neglecta_res1yr_20161207_1609/erd_v5_Sturnella_neglecta_predictors_test.txt'

target_col = os.path.basename(os.path.dirname(txt)).split('_res')[0]
out_png = '/vol/v2/stem/ebird/results/%s_obs_per_lc_2011_eval.png' % target_col
lc_path = '/vol/v1/proj/lst/outputs/models/randomforest/rfprediction_mosaic/yearly/lst_run1_prediction_voting_lulc_RF_mosaic_2011.bsq'
xlabels = ['water', 'snow', 'urban', 'bare', 'decid.', 'conif.', 'shrub']
obs_per_lc(txt, lc_path, target_col, out_png, xlabels)#'''


'''model_txt = '/vol/v2/stem/ebird/results/all_species_metrics.csv'
search_dir = '/vol/v2/stem/ebird/models'
constant_vars = ['EFFORT_HRS',
 'YEAR',
 'temp_max',
 'temp_min',
 'precip',
 'ecoregion',
 'population',
 'elevation']

species_list = ['c_guttatus', 'h_rustica', 'p_melanocephalus', 's_neglecta', 'z_leucophrys']

out_png = '/vol/v2/stem/ebird/results/importance/importance_per_res.png'
plot_importance(model_txt, species_list, search_dir, constant_vars, out_png)#'''
