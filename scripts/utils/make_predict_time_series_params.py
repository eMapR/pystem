# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 15:33:36 2017

@author: shooper
"""
import sys
import os
from stem import read_params

def main(predict_params, start_year, end_year, out_dir, txt_out_dir, n_jobs_pred=0, n_jobs_agg=0, confusion=False):
    
    param_dict, df_var_orig = read_params(predict_params)
    for k, v in param_dict.iteritems():
        param_dict[k] = v.replace('"','')
    
    param_basename = os.path.basename(predict_params)
    out_txt = os.path.join(txt_out_dir, param_basename)
    if not os.path.isdir(txt_out_dir):
        os.mkdir(txt_out_dir)
    
    for year in range(int(start_year), int(end_year) + 1):
        print 'Making params for year %s...' % year
        
        #################################################################################
        # jdb added 6/2/2017 - make a fresh copy of the original df_var otherwise line 32 won't work right
        df_var = df_var_orig.copy()     
        #################################################################################
        
        # Write the variable param table first
        band = year - 1983 # 1983 becuse gdal bands start at 1
        df_var.ix[df_var.data_band != 1, 'data_band'] = band
        df_var.data_band = df_var.data_band.astype(int)
        
        this_txt = out_txt.replace('.txt', '_%s.txt' % year)
        df_var.to_csv(this_txt, sep='\t')
        
        # Adjust a couple of variable values
        file_stamp = os.path.basename(param_dict['model_dir'].split('_')[0]) + '_' + str(year)
        param_dict['file_stamp'] = file_stamp
        param_dict['out_dir'] = os.path.join(out_dir, str(year))
        #if 'confusion_params' in param_dict and not confusion:
            #del param_dict['confusion_params']
        
        with open(this_txt, 'a') as txt:
            txt.write('\n')
            txt.write('model_dir; %s\n' % param_dict['model_dir'])
            txt.write('train_params; %s\n' % param_dict['train_params'])
            txt.write('mosaic_path; %s\n' % param_dict['mosaic_path'])
            txt.write('support_size; %s\n' % param_dict['support_size'])
            txt.write('n_tiles; %s\n' % param_dict['n_tiles'])
            txt.write('nodata; %s\n' % param_dict['nodata'])
            txt.write('out_dir; %s\n' % param_dict['out_dir'])
            txt.write('agg_stats; vote, pct_vote\n')
            txt.write('\nOptional Parameters\n')
            txt.write('file_stamp; %s\n' % param_dict['file_stamp'])
            if int(n_jobs_pred) != 0:
                n_jobs_pred = int(n_jobs_pred)
                txt.write('n_jobs_pred; %s\n' % n_jobs_pred)
            if int(n_jobs_agg) != 0:
                n_jobs_agg = int(n_jobs_agg)
                txt.write('n_jobs_agg; %s' % n_jobs_agg)
        
        print 'Params written to %s\n' % this_txt


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))