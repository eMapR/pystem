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
import warnings
import pandas as pd
import cPickle as pickle
from osgeo import gdal, ogr, gdal_array
from multiprocessing import Pool
import numpy as np

import stem_conus
import mosaic_by_tsa as mosaic
from lthacks import get_min_numpy_dtype, attributes_to_df
    
gdal.UseExceptions()


def parse_constant_vars(constant_vars):
    ''' helper function to isolate dict comprehension from exec statement'''
    
    var_dict = {k: int(v) for k, v in 
                [[i.replace(' ','') for i in item.replace(' ','').split(':')]
                for item in constant_vars.split(',')]}
    
    return var_dict

    
def main(params, inventory_txt=None, constant_vars=None, mosaic_shp=None, resolution=30, n_jobs_pred=0, n_jobs_agg=0, mosaic_nodata=0):
    
    inputs, df_var = stem_conus.read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])    
    df_var.data_band = [int(b) for b in df_var.data_band]#sometimes read as float

    try:
        support_size = [int(i) for i in support_size.split(',')]
        nodata = int(nodata)
        str_check = model_dir, mosaic_path, out_dir, train_params
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
    
    # Check that all the variables given were used in training and vice versa
    try:
        train_inputs, train_vars = stem_conus.read_params(train_params)
    except:
        raise NameError('train_params not specified or does not exist')
    train_vars = sorted(train_vars.index)
    pred_vars  = sorted(df_var.index)
    # Make sure vars are sorted alphabetically since they were for training
    df_var = df_var.reindex(pred_vars)
    
    unmatched_vars = [v for v in pred_vars if v not in train_vars]
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
    
    if mosaic_path.endswith('.shp'):
        mosaic_type = 'vector'
        if 'resolution' not in inputs:
            warnings.warn('Resolution not specified. Using default of 30...\n')
        mosaic_dataset = ogr.Open(mosaic_path)
        mosaic_ds = mosaic_dataset.GetLayer()
        min_x, max_x, min_y, max_y = mosaic_ds.GetExtent()
        xsize = int((max_x - min_x)/resolution)
        ysize = int((max_y - min_y)/resolution)
        prj = mosaic_ds.GetSpatialRef().ExportToWkt()
        x_res = resolution
        y_res = -resolution
        x_rot = 0
        y_rot = 0
        mosaic_tx, extent = stem_conus.tx_from_shp(mosaic_path, x_res, y_res)
        #df_tiles = attributes_to_df(mosaic_path)
    
    else:
        mosaic_type = 'raster'
        mosaic_ds = gdal.Open(mosaic_path)
        mosaic_tx = mosaic_ds.GetGeoTransform()
        xsize = mosaic_ds.RasterXSize
        ysize = mosaic_ds.RasterYSize
        prj = mosaic_ds.GetProjection()
        driver = mosaic_ds.GetDriver()
        m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    driver = gdal.GetDriverByName('gtiff')
    
    # If number of tiles not given, need to set it
    if 'n_tiles' not in inputs:
        print 'n_tiles not specified. Using default: 25 x 15 ...\n'
        n_tiles = 25, 15
    else:
        n_tiles = [int(i) for i in n_tiles.split(',')]
    df_tiles, df_tiles_rc, tile_size = stem_conus.get_tiles(n_tiles, xsize, ysize, mosaic_tx)
        
    predict_dir = os.path.join(out_dir, 'decisiontree_predictions')
    if not os.path.exists(predict_dir):
        os.mkdir(predict_dir)
    
    set_txt = glob.glob(os.path.join(model_dir, 'decisiontree_models/*support_sets.txt'))[0]
    df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    total_sets = len(df_sets)
    
    t0 = time.time()
    if 'n_jobs_pred' in inputs:
        n_jobs = int(n_jobs_pred)
        # Predict in parallel
        args = []
        t1 = time.time()
        print 'Predicting in parallel with %s jobs...' % n_jobs
        print 'Building args and making rasters of tile arrays...'
        for c, (set_id, row) in enumerate(df_sets.iterrows()):
            
            # Save rasters of tsa arrays ahead of time to avoid needing to pickle or fork mosaic
            coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
            '''if mosaic_type == 'vector':
                tsa_ar, tsa_off = mosaic.kernel_from_shp(mosaic_ds, coords, mosaic_tx, nodata=0)
            else:
                tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, 1, coords,
                                                        mosaic_tx, xsize, ysize,
                                                        nodata=nodata)
            set_mosaic_path = os.path.join(predict_dir, 'tsa_%s.tif' % set_id)
            tx_out = row.ul_x, mosaic_tx[1], mosaic_tx[2], row.ul_y, mosaic_tx[4], mosaic_tx[5]
            np_dtype = get_min_numpy_dtype(tsa_ar)
            gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(np_dtype)
            mosaic.array_to_raster(tsa_ar, tx_out, prj, driver, set_mosaic_path, gdal_dtype, silent=True)
            pct_progress = float(c + 1)/total_sets * 100
            sys.stdout.write('\rRetreived points for feature %s of %s (%%%.1f)' % (c + 1, total_sets, pct_progress))
            sys.stdout.flush()'''
            
            # Build list of args to pass to the Pool
            #tsa_off = stem_conus.calc_offset((mosaic_tx[0], mosaic_tx[3]), (tx_out[0], tx_out[3]), tx_out)
            args.append([coords, mosaic_type, mosaic_path, mosaic_tx, prj, nodata, c, total_sets, set_id, df_var, xsize, ysize, row.dt_file, nodata, np.uint8, constant_vars, predict_dir])
            #args.append([c, total_sets, set_id, df_var, set_mosaic_path, tsa_off, coords, 
                         #mosaic_tx, xsize, ysize, row.dt_file, nodata, np.uint8, 
                         #constant_vars, predict_dir])
        print '%.1f minutes\n' % ((time.time() - t1)/60)
        p = Pool(n_jobs)
        p.map(stem_conus.par_predict, args, 1)

    else:
        # Loop through each set and generate predictions
        for c, (set_id, row) in enumerate(df_sets.iterrows()):
            t1 = time.time()
            with open(row.dt_file, 'rb') as f: 
                dt_model = pickle.load(f)
            print '\nPredicting for set %s of %s' % (c + 1, total_sets)
            coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
            ar_predict = stem_conus.predict_set(set_id, df_var, mosaic_ds, coords, 
                                     mosaic_tx, xsize, ysize, dt_model, nodata,
                                     np.int16, constant_vars)        
            tx = coords.ul_x, x_res, x_rot, coords.ul_y, y_rot, y_res
            out_path = os.path.join(predict_dir, 'prediction_%s.tif' % set_id)
            np_dtype = get_min_numpy_dtype(ar_predict)
            gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(np_dtype)
            mosaic.array_to_raster(ar_predict, tx, prj, driver, out_path, gdal.GDT_Byte, nodata=nodata)
            print 'Total time for this set: %.1f minutes' % ((time.time() - t1)/60)
    
        #mosaic = None                  
    print '\nTotal time for predicting: %.1f hours\n' % ((time.time() - t0)/3600)#''' """
    
    #Aggregate predictions by tile and stitch them back together
    if not 'file_stamp' in inputs: file_stamp = os.path.basename(model_dir)
    
    t1 = time.time()
    agg_stats = [s.strip().lower() for s in agg_stats.split(',')]
    if 'n_jobs_agg' in inputs:
        n_jobs_agg = int(n_jobs_agg)
    
    if mosaic_type == 'vector':
        nodata_mask = mosaic_ds
    else:
        if 'mosaic_nodata' in inputs: mosaic_nodata = int(mosaic_nodata)
        nodata_mask = mosaic_ds.ReadAsArray() != mosaic_nodata
    
    pct_importance, df_sets = stem_conus.aggregate_predictions(n_tiles, ysize, xsize, nodata, nodata_mask, mosaic_tx, support_size, agg_stats, predict_dir, df_sets, out_dir, file_stamp, prj, driver, n_jobs_agg)
    #print 'Total aggregation time: %.1f hours\n' % ((time.time() - t0)/3600)
    mosaic_ds = None
    mosaic_dataset = None
    
    # Save the importance values
    importance = pd.DataFrame({'variable': pred_vars,
                               'pct_importance': pct_importance,
                               'index': range(len(pred_vars))
                               })
    importance.set_index('index', inplace=True)
    importance['rank'] = [int(r) for r in importance.pct_importance.rank(method='first', ascending=False)]
    out_txt = os.path.join(out_dir, '%s_importance.txt' % file_stamp)
    importance.to_csv(out_txt, sep='\t')#'''
    
    if 'confusion_params' in locals():
        import confusion_matrix as confusion

        ''' 
         Read the mean or vote back in '''
        if 'vote' in agg_stats:
            vote_path = os.path.join(out_dir, '%s_vote.tif')
            ar_vote = gdal.Open(vote_path)

        
            vote_dir = os.path.join(model_dir, 'evaluation_vote')
            mean_dir = os.path.join(model_dir, 'evaluation_mean')
            
            print '\nComputing confusion matrix for vote...'
            out_txt = os.path.join(vote_dir, 'confusion.txt')
            df_v = confusion.main(confusion_params, ar_vote, out_txt, match=True)
            vote_acc = df_v.ix['producer', 'user']
            vote_kap = df_v.ix['producer', 'kappa']
            '''try:
                out_txt = os.path.join(vote_dir, 'confusion_avg_kernel.txt')
                df_v_off = confusion.main(confusion_params, ar_vote, out_txt)
            except Exception as e:
                print e'''

                
        if 'mean' in agg_stats:
            mean_path = os.path.join(out_dir, '%s_mean.tif')
            ar_mean = gdal.Open(mean_path)
            print '\nGetting confusion matrix for mean...'
            out_txt = os.path.join(mean_dir, 'confusion.txt')
            df_m = confusion.main(confusion_params, ar_mean, out_txt, match=True)
            mean_acc = df_m.ix['user','producer']
            mean_kap = df_m.ix['user', 'kappa']
            '''try:
                out_txt = os.path.join(mean_dir, 'confusion_avg_kernel.txt')
                df_m_off = confusion.main(confusion_params, ar_mean, out_txt)
            except Exception as e:
                print e#'''


        if 'inventory_txt' in inputs:
            df_inv = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
            cols = ['vote_accuracy', 'vote_kappa']#, 'vote_mask', 'mean_accuracy', 'mean_kappa', 'vote_mask']
            df_inv.ix[file_stamp, cols] = vote_acc, vote_kap#, False, mean_acc, mean_kap, False
            df_inv.to_csv(inventory_txt, sep='\t')
        else:
            print '\n"inventory_txt" was not specified.' +\
            ' Model evaluation scores will not be recorded...'
            
        print ''
        if 'vote' in agg_stats:
            print 'Vote accuracy .............. ', vote_acc
            print 'Vote kappa ................. ', vote_kap
        if 'mean' in agg_stats:
            print 'Mean accuracy .............. ', mean_acc
            print 'Mean kappa ................. ', mean_kap
        
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

