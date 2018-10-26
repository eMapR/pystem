# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:32:12 2017

@author: shooper
"""
import os
import sys
import shutil
import pandas as pd

package_dir = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(package_dir)
from evaluation import confusion_matrix

def write_params(params_path, param_dict):
    
    lines = []
    lines.append('Parameters for confusion_matrix.py\n\n')
    lines.append('p_path; %s\n' % param_dict['p_path'])
    lines.append('t_path; %s\n' % param_dict['t_path'])
    lines.append('sample_txt; %s\n' % param_dict['sample_txt'])
    lines.append('target_col; %s\n' % param_dict['target_col'])
    lines.append('bins; %s\n' % param_dict['bins'])
    lines.append('p_nodata; %s\n' % param_dict['p_nodata'])
    lines.append('t_nodata; %s\n' % param_dict['t_nodata'])
    lines.append('\nOptional parameters\n')
    lines.append('out_txt; %s\n' % param_dict['out_txt'].replace('confusion.txt','revision_round3/confusion.txt'))
    lines.append('match; %s\n' % param_dict['match'])
    

    for i, l in enumerate(lines):
        if '/vol/v2/stem' in l and '/vol/v2/stem/caorwa' not in l:
            lines[i] = l.replace('/vol/v2/stem', '/vol/v2/stem/caorwa')
    #if param_dict['inventory_txt']:
    if param_dict['inventory_txt']:
        lines.append('inventory_txt; %s\n' % param_dict['inventory_txt'])
    lines.append('file_stamp; %s\n' % param_dict['file_stamp'])
        
    with open(params_path, 'w') as f:
        for l in lines: 
            f.write(l)
    
    print 'Params written to ', params_path, '\n'


def main(search_dir, out_dir, inventory_txt=None, t_path=None, agg_method='mean'):
    
    #stamps = [s for s in os.listdir(search_dir) if '.txt' not in s]
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    #stamps = pd.read_csv('/vol/v2/stem/caorwa/canopy/models/model_info_final.txt', sep='\t').stamp.tolist()
    stamps = pd.read_csv('/vol/v2/stem/caorwa/imperv/models/model_info_final.txt', sep='\t').stamp.tolist()
    
    for stamp in stamps:
        print stamp
        model_dir = os.path.abspath(os.path.join(search_dir, stamp))
        confusion_params = os.path.join(model_dir, 'confusion_matrix_params.txt')
        if not os.path.exists(confusion_params):
            print 'Skipping %s\n' % confusion_params
            continue
        
        param_dict = confusion_matrix.read_params(confusion_params)
        for k, v in param_dict.iteritems():
            param_dict[k] = v.replace('"','')
        param_dict['p_path'] = os.path.join(model_dir, stamp + '_%s.bsq' % agg_method)
        if t_path:
            param_dict['t_path'] = t_path

        '''sample_txt = param_dict['sample_txt']
        target_col = param_dict['target_col']
        bins = param_dict['bins']
        p_nodata = param_dict['p_nodata']
        t_nodata = param_dict['t_nodata']'''
        param_dict['match'] = 'best'
        param_dict['inventory_txt'] = inventory_txt
        param_dict['file_stamp'] = stamp
        
        eval_dir = os.path.join(model_dir, 'evaluation_%s' % agg_method)
        if not os.path.exists(eval_dir):
            print 'Eval dir does not exist for %s. Skipping...\n' % stamp
            continue
        out_txt = os.path.join(eval_dir, 'confusion.txt')
        param_dict['out_txt'] = out_txt
        
        params_path = os.path.join(out_dir, 'confusion_params_%s_%s.txt' % (stamp, agg_method))        
        write_params(params_path, param_dict)
        shutil.copy2(params_path, model_dir)
        
        '''this_srch_str = os.path.join(model_dir, 'train_stem_regressor_params.txt')
        matched = glob.glob(this_srch_str)
        if len(matched) == 0:
            print 'No training param file for ', stamp
            continue
        this_train_params = matched[0]
        t_param_dict, _ = stem.read_params(this_train_params)
        
        this_srch_str = os.path.join(model_dir, 'predict_stem_params.txt')
        matched = glob.glob(this_srch_str)
        if len(matched) == 0:
            print 'No predict param file for ', stamp
            continue
        this_predict_params = matched[0]
        p_param_dict, _ = stem.read_params(this_predict_params)'''
        

if __name__ == '__main__':
    
    sys.exit(main(*sys.argv[1:]))
        