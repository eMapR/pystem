import os
import sys
import shutil
import re
import fnmatch
import pandas as pd
import numpy as np
from datetime import datetime

import randomforest as forest

def main(params):

    # Read params and make variables from text
    inputs = forest.read_params(params)
    for i in inputs:
        #import pdb; pdb.set_trace()
        exec ("{0} = str({1})").format(i, inputs[i])

    # Check that variables were specified in params
    try:
        str_check = sample_txt, target_col, var_txt, out_dir
    except NameError as e:
        print ''
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)

    # Make optional numeric arguments numeric
    if 'n_trees' in locals():
        n_trees = int(n_trees)
    else:
        n_trees = 200
    if 'n_jobs' in locals():
        n_jobs = int(n_jobs)
    else:
        n_jobs = 1
    if 'max_depth' in locals():
        max_depth = int(max_depth)
    else:
        max_depth=None

    # Raise an error if var_txt doesn't exist. Otherwise, just read it in
    if not os.path.exists(var_txt):
        print ''
        msg = 'var_text path specified does not exist:\n%s\n\n' % var_txt
        raise IOError(msg)
    df_var = pd.read_csv(var_txt, sep='\t', index_col='var_name')
    
    # Make the output directory
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    if not 'out_dirname' in locals(): out_dirname = target_col
    stamp = '{0}_{1}_{2}'.format(out_dirname, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it already exists
    shutil.copy2(params, out_dir) #Copy the params so the parameters used are saved
    shutil.copy2(sample_txt, out_dir)

    # Read in training samples
    df_train = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')

    # Check that df_train has exactly the same columns as variables specified in df_vars
    train_columns = df_train.columns.tolist()
    unmatched_vars = [v for v in df_var.index if v not in train_columns]
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in sample_txt but specified in params:\n' + unmatched_str
        raise NameError(msg)

    # Sort the predictors in alphabetical order so that train columns can be in the same order as the predict array when
    #   predicting later on
    predict_cols = sorted(np.unique([c for c in df_train.columns if c in df_var.index]))
    predict_cols = [c for c in predict_cols if c in df_var.index]
    if target_col in predict_cols: predict_cols.remove(target_col)
    df_var = df_var.sort_index()
    if 'constant_vars' in inputs:
        constant_vars = sorted([i.strip() for i in constant_vars.split(',')])
        unmatched_vars = [v for v in constant_vars if v not in train_columns]
        if len(unmatched_vars) != 0:
            unmatched_str = '\n'.join(unmatched_vars)
            msg = 'Columns not in sample_txt but specified in params:\n' + unmatched_str
            raise NameError(msg)
        predict_cols += constant_vars

    x_train = df_train.reindex(columns=predict_cols)
    y_train = df_train[target_col]
    rf_model = forest.train_rf_regressor(x_train, y_train,
                                         ntrees=n_trees,
                                         njobs=n_jobs,
                                         max_depth=max_depth)
    if 'constant_vars' in inputs:
        for v in constant_vars:
            df_var = df_var.append(pd.Series(name=v))
    importance = rf_model.feature_importances_
    df_var['importance'] = importance
    df_var['rank'] = [int(r) for r in df_var.importance.rank(method='first', ascending=False)]
    out_txt = os.path.join(out_dir, '%s_importance.txt' % stamp)

    rf_path = os.path.join(out_dir, 'regressor_model_%s' % stamp)
    forest.save_rfmodel(rf_model, rf_path)
    oob_score = round(rf_model.oob_score_, 3)
    out_var_txt = os.path.join(out_dir, os.path.basename(var_txt))
    df_var.to_csv(out_var_txt, sep='\t')

    # Record params in inventory text file
    if 'inventory_txt' in inputs:
        df_inv = pd.read_csv(inventory_txt, sep='\t')
        cols = df_inv.columns
        try:
            res = int(re.search('[0-9]{1,2}', out_dirname).group())
        except:
            res = None
        df_inv = df_inv.append(pd.DataFrame([{'stamp': stamp,
                                              'temporal_res': res,
                                              'oob_score': oob_score, 
                                              'auc': None, 
                                              'rmse': None,
                                              'rmse_n': None,
                                              'rmse_p': None, 
                                              'n_samples': len(df_train),
                                              'n_trees': n_trees
                                              }]),
                               ignore_index=True)
        df_inv = df_inv.reindex(columns=cols)
        existing_models = fnmatch.filter(os.listdir(os.path.dirname(out_dir)), '*res*')
        df_inv = df_inv[df_inv.stamp.isin(existing_models)]
        df_inv.to_csv(inventory_txt, sep='\t', index=False)


    print 'Random Forest Regressor model written to:\n', rf_path
    print '\nOOB score: ', oob_score
    print 'Relative importance:'
    print df_var.importance.sort_values(ascending=False)


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))



