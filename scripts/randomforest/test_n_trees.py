import os
import sys
import shutil
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import style
from datetime import datetime

import randomforest as forest

import warnings
warnings.filterwarnings("ignore")


def test(out_dir, x_train, y_train, max_trees, step, min_trees=None, silent=False):
    ''' Make a plot for sample_txt of number of trees vs OOB error rates '''

    if not silent: print 'Testing OOB error rate per number of trees...'
    
    oob_errors = []
    if not min_trees:
        min_trees = step
    n_trees = range(min_trees, max_trees + 1, step)
    n_tests = len(n_trees)

    for i, n in enumerate(n_trees):
        rf_model = forest.train_rf_regressor(x_train, y_train, ntrees=n, njobs=1)
        oob_errors.append(1 - rf_model.oob_score_)
        if not silent:
            print 'Testing %s of %s models with %s trees: %.3f' % (i + 1, n_tests, n, rf_model.oob_score_)

    plt.plot(n_trees, oob_errors, '-')
    plt.axis([0, max_trees, 0, 1])
    plt.xlabel('Number of Trees')
    plt.ylabel('Out of Bag Error Rate')
    plt.title('Number of Decision Trees vs. OOB Error Rate')
    out_png = os.path.join(out_dir, 'ntrees_vs_oob_error.png')
    plt.savefig(out_png)
    plt.clf()

    if not silent: print 'Plot PNG written to : ', out_png, '\n'
    
    return [n_trees, oob_errors]


def main(params, constant_vars=[], silent=False, return_results=True):

    # Read params and make variables from text
    inputs = forest.read_params(params, silent=silent)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])

    # Check that variables were specified in params
    try:
        str_check = sample_txt, target_col, var_txt
        max_trees = int(max_trees)
        step = int(step)
    except NameError as e:
        print ''
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)

    # Try to read var_txt
    if not os.path.exists(var_txt):
        print ''
        msg = 'var_text path specified does not exist:\n%s\n\n' % var_txt
        raise IOError(msg)
    df_var = pd.read_csv(var_txt, sep='\t', index_col='var_name')
    
    if 'constant_vars' in inputs:
        constant_vars = sorted([i.strip() for i in constant_vars.split(',')])
        
    df_train = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    predict_cols = sorted(np.unique([c for c in df_train.columns for v in df_var.index if v in c] + constant_vars))
    
    y_train = df_train[target_col]
    x_train = df_train.reindex(columns=predict_cols)
    
    
    out_dir = os.path.dirname(sample_txt)
    scores = test(out_dir, x_train, y_train, max_trees, step, silent=silent)
    shutil.copy2(var_txt, out_dir)
    
    if return_results:
        return scores


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))



