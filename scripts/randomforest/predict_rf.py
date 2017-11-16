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
from multiprocessing import Pool

import randomforest as forest
from stem import evaluate_ebird, evaluate_by_lc
import stem
import mosaic_by_tsa as mosaic

matplotlib.style.use('ggplot')

def parse_constant_vars(constant_vars):
    ''' helper function to isolate dict comprehension from exec statement'''
    
    var_dict = {k: int(v) for k, v in 
                [[i.replace(' ','') for i in item.replace(' ','').split(':')]
                for item in constant_vars.split(',')]}
    
    return var_dict


def main(params, n_pieces=False, ydims=None, constant_vars=None, year='', agg_method=None):

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
    '''if 'constant_vars' in inputs:
        constant_vars = parse_constant_vars(constant_vars)
        #year = constant_vars['YEAR']
        year = 2012
        pred_constants = sorted(constant_vars.keys())
    else:
        df_var.search_str = [s.format(2007) for s in df_var.search_str]'''

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
    if 'agg_method' in inputs:
        agg_method = inputs['agg_method']
        
    # Get mask and raster info
    ds = gdal.Open(mask_path)
    ar = ds.ReadAsArray()
    nodata_mask = ar != 0
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    tx = ds.GetGeoTransform()
    prj = ds.GetProjection()
    driver = gdal.GetDriverByName('gtiff')
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    

    # Predict
    #print 'Predicting with %s processors... %s' % (rf_model.n_jobs, time.ctime(time.time()))
    t1 = time.time()
    predict_pieces = []

    
    if 'n_tiles' not in inputs:
        print 'n_tiles not specified. Using default: 25 x 15 ...\n'
        n_tiles = 25, 15
    else:
        n_tiles = [int(i) for i in n_tiles.split(',')]
        
    if 'n_tiles' in inputs:
        df_tiles, df_tiles_rc, tile_size = stem.get_tiles(n_tiles, xsize, ysize, tx)
        empty_tiles = []
        ar_out = np.full((ysize, xsize), nodata, dtype=np.uint8)
        tile_dir = os.path.join(out_dir, 'predict_tiles')
        if not os.path.isdir(tile_dir):
            os.mkdir(tile_dir)
        for i, (ind, tile_coords) in enumerate(df_tiles.iterrows()):
            print 'Predicting for tile %s of %s...' % (i + 1, len(df_tiles))
            t1 = time.time()
            coords = tile_coords[['ul_x', 'ul_y', 'lr_x', 'lr_y']].tolist()
            tsa_ar, tsa_off = mosaic.extract_kernel(ds, 1, coords, tx, xsize, ysize, nodata=nodata)
            tsa_mask = tsa_ar == 0
            if tsa_mask.all():
                print 'Tile %s empty. Skipping...' % ind
                continue
            tsa_ar[tsa_mask] = nodata
            # Get the ids of TSAs this kernel covers
            tsa_ids = np.unique(tsa_ar)
            #tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=nodata]
            tsa_strs = [str(tsa) for tsa in tsa_ids if tsa!=nodata]
            array_shape = tsa_ar.shape
        
            # Get an array of predictors where each column is a flattened 2D array of a
            #   single predictor variable
            temp_nodata = -9999
            ar_predictors = stem.get_predictors(df_var, tx, tsa_strs, tsa_ar, coords, tsa_mask, temp_nodata, 1)
            nodata_mask = ~ np.any(ar_predictors==temp_nodata, axis=1)
            predictors = ar_predictors[nodata_mask]
            t2 = time.time()
            if agg_method == 'mode':
                args = []
                for dt in rf_model.estimators_:
                    args.append([dt, predictors])
                pool = Pool(rf_model.n_jobs)
                t3 = time.time()
                dt_predictions = np.vstack(pool.map(forest.par_predict_from_dt, args, 1))
                print 'Prediction time: %.1f minutes' % ((time.time() - t3)/60)
                t3 = time.time()
                predictions = stem.mode(dt_predictions, axis=0)
                print 'Aggregation time:  %.1f minutes' % ((time.time() - t3)/60)
                del dt_predictions
                t3 = time.time()
                pool.close()
                pool.join()
                print 'Closing time:  %.1f minutes' % ((time.time() - t3)/60)
            else:
                predictions = rf_model.predict(ar_predictors[nodata_mask])
            print 'Prediction time: %.1f minutes' % ((time.time() - t2)/60)
            
            ar_tile = np.full(ar_predictors.shape[0], nodata, dtype=np.uint8)
            ar_tile[nodata_mask] = predictions.astype(np.uint8)
            ul_r, lr_r, ul_c, lr_c = df_tiles_rc.ix[ind]
            ar_out[ul_r : lr_r, ul_c : lr_c] = ar_tile.reshape(array_shape)
            tx_tile = tile_coords.ul_x, x_res, x_rot, tile_coords.ul_y, y_rot, y_res
            mosaic.array_to_raster(ar_tile.reshape(array_shape), tx_tile, prj, driver, os.path.join(tile_dir, 'tile_%s.tif' % ind), dtype=gdal.GDT_Byte, nodata=nodata)
            print 'Total time for this piece: %.1f minutes\n' % ((time.time() - t1)/60)
            #del ar_predictors, nodata_mask, ar_prediction'''
        #ar_prediction = np.concatenate(predict_pieces)
        #del predict_pieces
        '''ar_out = np.full((ysize, xsize), nodata, dtype=np.uint8)
        for ind, tile_coords in df_tiles_rc.iterrows():
            if ind in empty_tiles:
                continue
            ul_r, lr_r, ul_c, lr_c = tile_coords
            tile_file = os.path.join(tile_dir, 'tile_%s.tif' % ind)
            if not os.path.exists(tile_file):
                continue
            ds_t = gdal.Open(tile_file)
            ar_tile = ds_t.ReadAsArray()
            t_ulx = df_tiles.ix[ind, ['ul_x', 'ul_y']]
            ar_out[ul_r : lr_r, ul_c : lr_c] = ar_tile'''
        
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
    out_path = os.path.join(out_dir, '%s_rf_vote.tif' % stamp)
    #ar_prediction = ar_prediction.reshape(ysize, xsize)
    if constant_vars: 
        out_path = out_path.replace('.tif', '_yr%s.tif' % year )
    forest.array_to_raster(ar_out, tx, prj, driver, out_path, gdal.GDT_Byte, nodata)#"""
    # Delete the tiles
    shutil.rmtree(tile_dir)
    ds = None
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

    print '\nTotal runtime: %.1f minutes' % ((time.time() - t0)/60)
    
    return out_path


if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))




