# -*- coding: utf-8 -*-
"""
Predict from a spatiotemporal exploratory model 

@author: Sam Hooper, samhooperstudio@gmail.com

"""
import glob
import time
import os
import sys
import shutil
import pandas as pd
import cPickle as pickle
from osgeo import gdal
from multiprocessing import Pool
import numpy as np

import stem
import mosaic_by_tsa as mosaic
    
gdal.UseExceptions()


def parse_constant_vars(constant_vars):
    ''' helper function to isolate dict comprehension from exec statement'''
    
    var_dict = {k: int(v) for k, v in 
                [[i.replace(' ','') for i in item.replace(' ','').split(':')]
                for item in constant_vars.split(',')]}
    
    return var_dict

    
def main(params, inventory_txt=None, constant_vars=None):
    
    inputs, df_var = stem.read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])    
    df_var.data_band = [int(b) for b in df_var.data_band]#sometimes read as float

    try:
        n_tiles = [int(i) for i in n_tiles.split(',')]
        support_size = [int(i) for i in support_size.split(',')]
        nodata = int(nodata)
        str_check = model_dir, mosaic_path, out_dir, train_params
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    # Check that all the variables given were used in training and vice versa
    try:
        train_inputs, train_vars = stem.read_params(train_params)
    except:
        raise NameError('train_params not specified or does not exist')
    train_vars = sorted(train_vars.index)
    pred_vars  = sorted(df_var.index)
    # Make sure vars are sorted alphabetically since they were for training
    df_var = df_var.reindex(pred_vars)
    # If constants were given, make a dict and make sure they match the training
    #  constants
    if 'constant_vars' in inputs:
        constant_vars = parse_constant_vars(constant_vars)
        pred_constants = sorted(constant_vars.keys())
        train_constants = [i.replace(' ', '') for i in train_inputs['constant_vars'].strip('"').split(',')]
        train_constants = sorted(train_constants)
    unmatched_vars = [v for v in pred_vars if v not in train_vars]
    if 'constant_vars' in inputs: 
        unmatched_vars += [v for v in pred_constants if v not in train_constants]
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in train params but specified in predict params:\n' + unmatched_str
        raise NameError(msg)
    unmatched_vars = [v for v in train_vars if v not in pred_vars]
    if 'constant_vars' in inputs:
        unmatched_vars += [v for v in train_constants if v not in pred_constants]
        pred_vars += pred_constants # Add here because it would screw with stuff upstream
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in predict params but specified in train params:\n' + unmatched_str
        raise NameError(msg)
    
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    else: print ('WARNING: out_dir already exists:\n%s\nAny existing files ' + \
    'will be overwritten...\n') % out_dir
    if not os.path.exists(os.path.join(out_dir, os.path.basename(params))):
        shutil.copy2(params, out_dir) #Copy the params for reference
    
    if 'confusion_params' in inputs: 
        conf_bn = os.path.basename(confusion_params)
        new_conf_path = os.path.join(out_dir, conf_bn)
        if not os.path.exists(new_conf_path):
            shutil.copy2(confusion_params, out_dir)
        confusion_params = new_conf_path
    
    if not os.path.exists(model_dir):
        sys.exit('model_dir does not exist:\n%s' % model_dir)
    if not os.path.exists(mosaic_path):
        sys.exit('mosaic_path does not exist:\n%s' % mosaic_path)
    
    mosaic_ds = gdal.Open(mosaic_path)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    xsize = mosaic_ds.RasterXSize
    ysize = mosaic_ds.RasterYSize
    prj = mosaic_ds.GetProjection()
    driver = mosaic_ds.GetDriver()
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    
    predict_dir = os.path.join(out_dir, 'decisiontree_predictions')
    if not os.path.exists(predict_dir):
        os.mkdir(predict_dir)
    
    set_txt = glob.glob(os.path.join(model_dir, 'decisiontree_models/*support_sets.txt'))[0]
    df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    total_sets = len(df_sets)
    
    t0 = time.time()
    if 'n_jobs' in inputs:

        # Predict in parallel
        n_jobs = int(n_jobs)
        args = []
        t1 = time.time()
        print 'Predicting in parallel with %s jobs...' % n_jobs
        print 'Building args and making rasters of TSA arrays...'
        for c, (set_id, row) in enumerate(df_sets.iterrows()):
            
            # Save rasters of tsa arrays ahead of time to avoid needing to pickle or fork mosaic_ds
            coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
            tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, 1, coords, mosaic_tx,
                            xsize, ysize, nodata=nodata)
            tsa_raster = os.path.join(predict_dir, 'tsa_%s.bsq' % set_id)
            tx_out = row.ul_x, mosaic_tx[1], mosaic_tx[2], row.ul_y, mosaic_tx[4], mosaic_tx[5]
            dtype_code = mosaic_ds.GetRasterBand(1).DataType
            mosaic.array_to_raster(tsa_ar, tx_out, prj, driver, tsa_raster, stem.get_gdal_dtype(dtype_code), silent=True)
            
            # Build list of args to pass to the Pool
            tsa_raster = os.path.join(predict_dir, 'tsa_%s.bsq' % set_id)
            ds = gdal.Open(tsa_raster)
            tsa_tx = ds.GetGeoTransform()
            ds = None
            tsa_off = stem.calc_offset((mosaic_tx[0], mosaic_tx[3]), (tsa_tx[0], tsa_tx[3]), tsa_tx)
            args.append([c, total_sets, set_id, df_var, tsa_raster, tsa_off, coords, 
                         mosaic_tx, xsize, ysize, row.dt_file, nodata, np.uint8, 
                         constant_vars, predict_dir])
        print '%.1f minutes\n' % ((time.time() - t1)/60)
        p = Pool(n_jobs)
        p.map(stem.par_predict, args, 1)
            
    
    else:
        # Loop through each set and generate predictions
        for c, (set_id, row) in enumerate(df_sets.ix[1043:].iterrows()):
            t1 = time.time()
            with open(row.dt_file, 'rb') as f: 
                dt_model = pickle.load(f)
            print '\nPredicting for set %s of %s' % (c + 1, total_sets)
            coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
            ar_predict = stem.predict_set(set_id, df_var, mosaic_ds, coords, 
                                     mosaic_tx, xsize, ysize, dt_model, nodata,
                                     np.int16, constant_vars)        
            tx = coords.ul_x, x_res, x_rot, coords.ul_y, y_rot, y_res
            out_path = os.path.join(predict_dir, 'prediction_%s.bsq' % set_id)
            mosaic.array_to_raster(ar_predict, tx, prj, driver, out_path, gdal.GDT_Byte, nodata=nodata)
            print 'Total time for this set: %.1f minutes' % ((time.time() - t1)/60)
    
        #mosaic_ds = None                  
    print '\nTotal time for predicting: %.1f hours\n' % ((time.time() - t0)/3600)#'''
    
    #Aggregate predictions by tile and stitch them back together
    if not 'file_stamp' in inputs: file_stamp = os.path.basename(model_dir)
    ar_vote, pct_importance, df_sets = stem.aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_ds, support_size, predict_dir, df_sets, out_dir, file_stamp, prj, driver, 0)
    #df_sets.to_csv(set_txt, sep='\t')'''
    mosaic_ds = None
    
    '''# Save the importance values
    importance = pd.DataFrame({'variable': pred_vars,
                               'pct_importance': pct_importance,
                               'index': range(len(pred_vars))
                               })
    importance.set_index('index', inplace=True)
    importance['rank'] = [int(r) for r in importance.pct_importance.rank(method='first', ascending=False)]
    out_txt = os.path.join(out_dir, '%s_importance.txt' % file_stamp)
    importance.to_csv(out_txt, sep='\t')#'''
    
    '''ds = gdal.Open(os.path.join(model_dir, '%s_vote.bsq' % file_stamp))
    ar_vote = ds.ReadAsArray()
    ds = None
    ds = gdal.Open(os.path.join(model_dir, '%s_mean.bsq' % file_stamp))
    ar_mean = ds.ReadAsArray()
    ds = None#'''
    
    if 'confusion_params' in locals():
        import confusion_matrix as confusion
        
        vote_dir = os.path.join(model_dir, 'evaluation_vote')
        mean_dir = os.path.join(model_dir, 'evaluation_mean')
        
        '''print '\nComputing confusion matrix for vote...'
        out_txt = os.path.join(vote_dir, 'confusion.txt')
        print confusion_params
        df_v = confusion.main(confusion_params, ar_vote, out_txt, match=True)
        try:
            out_txt = os.path.join(vote_dir, 'confusion_avg_kernel.txt')
            df_v_off = confusion.main(confusion_params, ar_vote, out_txt)
        except Exception as e:
            print e'''
        
        print '\nGetting confusion matrix for mean...'
        out_txt = os.path.join(mean_dir, 'confusion.txt')
        df_m = confusion.main(confusion_params, ar_mean, out_txt, match=True)
        '''try:
            out_txt = os.path.join(mean_dir, 'confusion_avg_kernel.txt')
            df_m_off = confusion.main(confusion_params, ar_mean, out_txt)
        except Exception as e:
            print e#'''
        
        vote_acc = df_v.ix['producer', 'user']
        vote_kap = df_v.ix['producer', 'kappa']
        #mean_acc = df_m.ix['user','producer']
        #mean_kap = df_m.ix['user', 'kappa']

        if 'inventory_txt' in inputs:
            df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
            cols = ['vote_accuracy', 'vote_kappa']#, 'vote_mask', 'mean_accuracy', 'mean_kappa', 'vote_mask']
            df_inv.ix[file_stamp, cols] = vote_acc, vote_kap#, False, mean_acc, mean_kap, False
            df_inv.to_csv(inventory_txt, sep='\t')
        else:
            print '\n"inventory_txt" was not specified.' +\
            ' Model evaluation scores will not be recorded...'
            
        print ''
        print 'Vote accuracy .............. ', vote_acc
        print 'Vote kappa ................. ', vote_kap
        #print 'Mean accuracy .............. ', mean_acc
        #print 'Mean kappa ................. ', mean_kap
        
    else:
        print '\n"confusion_params" was not specified.' +\
            ' This model will not be evaluated...' #'''
    
    print '\nTotal prediction runtime: %.1f\n' % ((time.time() - t0)/60)

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))#'''


''' ############# Testing ################ '''
#sample_txt = '/vol/v2/stem/canopy/samples/canopy_sample3000_20160122_1600_predictors.txt'
#target_col = 'value'
#mosaic_path = '/vol/v1/general_files/datasets/spatial_data/CAORWA_TSA_lt_only.bsq'
#tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'
#cell_size = (300000, 200000)
#support_size = (400000, 300000)
#sets_per_cell = 10
#min_obs = 25
#pct_train = .63
#target_col = 'canopy'
#n_tiles = 10
#out_dir = '/vol/v2/stem/canopy/models/'
#set_id, ar, df_sets, df_train = main(sample_txt, target_col, mosaic_path, cell_size, support_size, sets_per_cell, min_obs, pct_train, target_col, n_tiles, out_dir)

#params = '/vol/v2/stem/param_files/build_stem_params_nomse.txt'
#predictions, df_sets, df_train = main(params)
#stuff = main(params)
'''set_txt = '/vol/v2/stem/canopy/outputs/canopy_20160212_2016/decisiontree_models/canopy_20160212_2016_support_sets.txt'
df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
tsa_ar = predict_set_from_disk(df_sets, 341, params)'''

'''tile_size = [size * 30 for size in tile_size]
shp = '/vol/v2/stem/extent_shp/orwaca.shp'
coords, extent = gsrd.get_coords(shp)
gsrd.plot_sets_on_shp(coords, 900, sets[20][1], (400000, 300000), df_tiles.ix[sets[20][0]], tile_size)'''
#for s in sets:
#    print 'Number of sets: ', len(s)
#    coords, extent = gsrd.get_coords(shp)
#    gsrd.plot_sets_on_shp(coords, 500, s, support_size)
