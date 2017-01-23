import glob
import time
import os
import sys
import shutil
import pandas as pd
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib
from osgeo import gdal

import randomforest as forest
from stem import evaluate_ebird, evaluate_by_lc

matplotlib.style.use('ggplot')

def parse_constant_vars(constant_vars):
    ''' helper function to isolate dict comprehension from exec statement'''
    
    var_dict = {k: int(v) for k, v in 
                [[i.replace(' ','') for i in item.replace(' ','').split(':')]
                for item in constant_vars.split(',')]}
    
    return var_dict


def main(params, n_pieces=False, ydims=None, constant_vars=None, year=''):

    t0 = time.time()
    print 'Predicting Random Forest... %s\n' % time.ctime(t0)

    # Set optional params to default:
    split_predictors = False

    # Read params and make variables from text
    inputs = forest.read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])

    # Check that variables were specified in params
    try:
        nodata = int(nodata)
        str_check = train_params, rf_path, mask_path, out_dir
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)

    # Raise an error if the var_txt path doesn't exist. Otherwise, just read it in
    train_dict = forest.read_params(train_params)
    train_txt_bn = os.path.basename(train_dict['var_txt'][:-1])
    if 'var_txt' not in locals():
        var_txt = os.path.join(os.path.dirname(rf_path), train_txt_bn)
    if not os.path.exists(var_txt):
        print ''
        msg = 'Could not find var_txt:\n%s\n' % var_txt
        raise IOError(msg)
    df_var = pd.read_csv(var_txt, sep='\t', index_col='var_name')

    # Make sure vars are sorted alphabetically since they were for training
    pred_vars  = sorted(df_var.index)
    df_var = df_var.reindex(pred_vars)
    if 'constant_vars' in inputs:
        constant_vars = parse_constant_vars(constant_vars)
        #year = constant_vars['YEAR']
        year = 2012
        pred_constants = sorted(constant_vars.keys())
    else:
        df_var.search_str = [s.format(2007) for s in df_var.search_str]

    #out_dir = os.path.dirname(out_raster)
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    else: print ('WARNING: out_dir already exists:\n%s\nAny existing files ' + \
    'will be overwritten...\n') % out_dir
    new_params = os.path.join(out_dir, os.path.basename(params))
    shutil.copy2(params, new_params.replace('.txt', '_%s.txt' % year))

    # Load the Random Forest model
    print 'Loading the RandomForest model from \n%s... \n%s\n' % (rf_path, time.ctime(time.time()))
    if not os.path.exists(rf_path):
        raise IOError('%s does not exist' % rf_path)
    with open(rf_path) as f:
        rf_model = pickle.load(f)
    n_features = rf_model.n_features_
    n_vars = len(df_var.index.tolist())
    if 'constant_vars' in inputs: 
        n_vars += len(pred_constants)
    if n_features != n_vars:
        print df_var.index.tolist() + pred_constants
        sys.exit(('\nKeyError: Number of features of the random forest model does not match the number of variables in df_var.' +\
            '\nNumber of features of the model: {0} \nNumber of variables in var_txt: {1}' + \
            '\nCheck that all predictors for used in var_txt to train the model are in this var_txt ' +\
            '\nPath of Random Forest model: {2}\nPath of var_txt: {3}').format(n_features, n_vars, rf_path, var_txt))
        #"""

    # Get mask and raster info
    ds = gdal.Open(mask_path)
    ar = ds.ReadAsArray()
    mask = ar != 0
    ar = None
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    driver = gdal.GetDriverByName('gtiff')
    ul_x, x_res, _, ul_y, _, y_res = tx
    ds = None

    # Predict
    'Predicting with %s processors... %s' % (rf_model.n_jobs, time.ctime(time.time()))
    t1 = time.time()
    predict_pieces = []
    
    if n_pieces:
        # assumes predictors all line up and have same dimensions'''
        if 'mask_path' not in inputs: 
            raise NameError('mask_path not specified')
        # Figure out the y dimension of each piece
        n_pieces = int(n_pieces)
        piece_ysize = ysize/n_pieces
        upper_ydim = range(0, ysize, piece_ysize)
        lower_ydim = range(piece_ysize, ysize, piece_ysize)
        lower_ydim[-1] = ysize
        ydims = zip(upper_ydim, lower_ydim)
        for i, yd in enumerate(ydims):
            print 'Predicting for piece %s of %s...' % (i + 1, n_pieces)
            t1 = time.time()
            ar_predictors, nodata_mask = forest.get_predictors(df_var, nodata, yd, constant_vars)
            t2 = time.time()
            predictions = rf_model.predict(ar_predictors)
            print 'Prediction time: %.1f minutes' % ((time.time() - t2)/60)
            ar_prediction = np.full(nodata_mask.shape[0], nodata, dtype=np.uint8)
            ar_prediction[nodata_mask] = (predictions * 100).astype(np.uint8)
            predict_pieces.append(ar_prediction)
            print 'Total time for this piece: %.1f minutes\n' % ((time.time() - t1)/60)
            del ar_predictors, nodata_mask, ar_prediction
        ar_prediction = np.concatenate(predict_pieces)
        del predict_pieces
    
    else:
        ar_predictors, nodata_mask = forest.get_predictors(df_var, nodata)
        # If the predictions are too large (i.e. cause memory errors), split the predictor array into pieces and predict
        #   separately, then stack them back together
        if split_predictors:
            split_predictors = int(split_predictors)
            predictions = []
            for i, p in enumerate(np.array_split(ar_predictors, split_predictors)):
                t1 = time.time()
                print '\nPredicting for %s of %s pieces of the final array...' % (i + 1, split_predictors)
                predictions.append(rf_model.predict(p))
                print '%.1f minutes' % ((time.time() - t1)/60)
            predictions = np.concatenate(predictions)
            print ''
        else:
            print 'Predicting in one chunk...'
            predictions = rf_model.predict(ar_predictors)
        ar_prediction = np.full(nodata_mask.shape[0], nodata, dtype=np.float32)
        ar_prediction[nodata_mask] = predictions
        del ar_predictors, predictions

    # Save the prediction array to disk
    stamp = os.path.basename(out_dir)
    out_path = os.path.join(out_dir, 'final_%s.tif' % stamp)
    ar_prediction = ar_prediction.reshape(ysize, xsize)
    if constant_vars: 
        out_path = out_path.replace('.tif', '_yr%s.tif' % year )
    forest.array_to_raster(ar_prediction, tx, prj, driver, out_path, gdal.GDT_Byte, nodata)#"""
    
    '''stamp = os.path.basename(out_dir)
    path = os.path.join(out_dir, 'final_%s_yr2011.tif' % stamp) 
    stamp = os.path.basename(os.path.dirname(path))
    ds = gdal.Open(path)
    ar_prediction = ds.ReadAsArray()
    ds = None#'''
    

    if 'test_params' in inputs:
        #df_test = pd.read_csv(test_samples, sep='\t', index_col='obs_id')
        print '\nEvaluating the model...'
        t1 = time.time()
        test_dict = forest.read_params(test_params)
        for i in test_dict:
            exec ("{0} = str({1})").format(i, test_dict[i])
            
        if 'n_trials' in test_dict: 
            n_trials = int(n_trials)
        else:
            'n_trials not specified. Setting default to 50...\n'
            n_trials = 50
        if 'year' in test_dict: 
            year = int(year)
        else:
            year = None
        cell_size = [int(i) for i in cell_size.split(',')]
        n_per_cell = int(n_per_cell)
        param_bn = os.path.basename(test_params)
        shutil.copy2(test_params, 
                     os.path.join(out_dir, 
                                  param_bn.replace('.txt', '_%s.txt' % year))
                    )
        
        df, samples, roc_curves = evaluate_ebird(sample_txt, ar_prediction, tx,
                                                 cell_size, target_col, n_per_cell,
                                                 n_trials, year)
                                                 
        if 'inventory_txt' in test_dict:
            score_cols = sorted(df.columns)
            df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp') 
            for col in score_cols:
                score_mean = df[col].mean()
                df_inv.ix[stamp, col] = score_mean
                print 'Average %s: %2.3f' % (col.upper(), score_mean) 
            df_inv.to_csv(inventory_txt, sep='\t')
        out_txt = os.path.join(out_dir, '{0}_eval_{1}.txt'.format(stamp, year))
        df.to_csv(out_txt, sep='\t', index=False)
        samples.to_csv(out_txt.replace('.txt', '_samples.txt'), sep='\t')
        print '\nTotal eval time: %.1f minutes\n' % ((time.time() - t1)/60)
    else:
        print '\nEither "test_samples" or "inventory_txt" was not specified.' +\
            ' This model will not be evaluated...'
    
    if len(roc_curves) > 0:
        for fpr, tpr, thresholds in roc_curves:
            plt.plot(fpr, tpr, 'k', alpha=.1)
        out_png = os.path.join(out_dir, '{0}_roc_curve_{1}.png'.format(stamp, year))
        plt.savefig(out_png)
    
    if 'lc_path' in test_dict:
        '''df_lc = evaluate_by_lc(samples, ar_prediction, lc_path, target_col)
        out_txt = os.path.join('/vol/v2/stem/ebird/results/performance_by_lc', '{0}_eval_{1}_land_cover.txt'.format(stamp, year))
        df_lc.to_csv(out_txt, sep='\t')'''
        
        #df_samples = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
        df_lc = evaluate_by_lc(samples, ar_prediction, lc_path, target_col)
        out_txt = os.path.join(out_dir, '{0}_eval_{1}_land_cover_all_samples.txt'.format(stamp, year))
        df_lc.to_csv(out_txt, sep='\t')

    print '\nTotal runtime: %.1f minutes' % ((time.time() - t0)/60)


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))




