# -*- coding: utf-8 -*-
"""
Predict from a spatiotemporal exploratory model 

@author: Sam Hooper, samhooperstudio@gmail.com

"""
import os
import sys
import glob
import re
import time
import fnmatch
import subprocess
import shutil
import warnings
import sqlalchemy
import pandas as pd
import cPickle as pickle
from osgeo import gdal, ogr, gdal_array
from sklearn.externals import joblib
from multiprocessing import Pool
import numpy as np

import stem
from evaluation.evaluation import feature_to_mask
import mosaic_by_tsa as mosaic
    
gdal.UseExceptions()

'''
TODO: 
-add support for generic on-the-fly tiles if user doesn't want to use 
mosaic_shp tiles
'''


def parse_constant_vars(constant_vars):
    ''' helper function to isolate dict comprehension from exec statement'''
    
    var_dict = {k: int(v) for k, v in 
                [[i.replace(' ','') for i in item.replace(' ','').split(':')]
                for item in constant_vars.split(',')]}
    
    return var_dict


def parallel_helper(n_jobs, df_var, tile_str, tile_ul, mosaic_tx, overlapping_sets, agg_stats, inds):
    return joblib.Parallel(n_jobs=n_jobs, verbose=5, backend="threading")(
joblib.delayed(stem.predict_pixel)
            (df_var, tile_str, tile_ul, mosaic_tx, pixel_inds, overlapping_sets, agg_stats) for pixel_inds in inds)

    

def main(params, inventory_txt=None, constant_vars=None, mosaic_shp=None, resolution=30, n_jobs=0, n_jobs_agg=0, mosaic_nodata=0, snap_coord=None, overwrite_tiles=False, tile_id_field='name'):
    inputs, df_var = stem.read_params(params)
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
        train_inputs, train_vars = stem.read_params(train_params)
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
    
    if not 'file_stamp' in inputs: file_stamp = os.path.basename(model_dir)
    db_path = os.path.join(model_dir, os.path.basename(model_dir) + '.db')
    if os.path.exists(db_path):
        engine = sqlalchemy.create_engine('sqlite:///%s' % db_path)
        with engine.connect() as con, con.begin():
            df_sets = pd.read_sql_table('support_sets', con, index_col='set_id')#'''
    else:
        set_txt = glob.glob(os.path.join(model_dir, 'decisiontree_models/*support_sets.txt'))[0]
        if not os.path.isfile(set_txt):
            raise IOError('No database or support set txt file found')
        df_sets = pd.read_csv(set_txt, sep='\t', index_col='set_id')
    
    if mosaic_path.endswith('.shp'):
        mosaic_type = 'vector'
        # if subset specified, clip the mosaic and set mosaic path to clipped shp
        if 'subset_shp' in inputs:
            out_shp_bn = os.path.basename(mosaic_path).replace('.shp', '_clipped.shp')
            out_shp = os.path.join(out_dir, out_shp_bn)
            cmd = 'ogr2ogr -clipsrc {clip_shp} {out_shp} {in_shp}'.format(clip_shp=subset_shp, out_shp=out_shp, in_shp=mosaic_path)
            subprocess.call(cmd, shell=True)#'''
            mosaic_path = out_shp
        mosaic_dataset = ogr.Open(mosaic_path)
        mosaic_ds = mosaic_dataset.GetLayer()
        min_x, max_x, min_y, max_y = mosaic_ds.GetExtent()
        if 'resolution' not in inputs:
            warnings.warn('Resolution not specified. Using default of 30...\n')
        # If subset specified, just get sets that overlap the subset
        if 'subset_shp' in inputs:
            mosaic_geom = ogr.Geometry(ogr.wkbMultiPolygon)
            for feature in mosaic_ds:
                mosaic_geom.AddGeometry(feature.GetGeometryRef())
            df_sets = stem.get_overlapping_sets(df_sets, mosaic_geom)
        xsize = int((max_x - min_x)/resolution)
        ysize = int((max_y - min_y)/resolution)
        prj = mosaic_ds.GetSpatialRef().ExportToWkt()
        x_res = resolution
        y_res = -resolution
        x_rot = 0
        y_rot = 0
        if 'snap_coord' in train_inputs:
            snap_coord = train_inputs['snap_coord'].replace('"','')
            snap_coord = [float(c) for c in snap_coord.split(',')]#'''
        mosaic_tx, extent = stem.tx_from_shp(mosaic_path, x_res, y_res, snap_coord=snap_coord)
        tiles = stem.attributes_to_df(mosaic_path) # Change to accept arbittary geometry
        
    else:
        mosaic_type = 'raster'
        mosaic_ds = gdal.Open(mosaic_path)
        mosaic_tx = mosaic_ds.GetGeoTransform()
        xsize = mosaic_ds.RasterXSize
        ysize = mosaic_ds.RasterYSize
        prj = mosaic_ds.GetProjection()
        driver = mosaic_ds.GetDriver()
        m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    #driver = gdal.GetDriverByName('gtiff')
        
    # If number of tiles not given, need to set it
    if 'n_tiles' not in inputs:
        print 'n_tiles not specified. Using default: 25 x 15 ...\n'
        n_tiles = 90, 40
    else:
        n_tiles = [int(i) for i in n_tiles.split(',')]
    #df_tiles, df_tiles_rc, tile_size = stem.get_tiles(n_tiles, xsize, ysize, mosaic_tx)
    
    total_sets = len(df_sets)
    t0 = time.time()
    last_dts = pd.Series()
    agg_stats = [s.strip().lower() for s in agg_stats.split(',')]
    n_jobs = int(n_jobs)
    tile_dir = os.path.join(out_dir, '_temp_tiles')
    #tile_dir = '/home/server/pi/homes/shooper/delete_test'
    if not os.path.isdir(tile_dir):
        os.mkdir(tile_dir)
    tile_path_template = os.path.join(tile_dir, 'tile_{tile_id}_%(stat)s.tif')
    n_tiles = len(tiles)
    
    if not overwrite_tiles:
        files = os.listdir(tile_dir)
        tile_files = pd.DataFrame(columns=agg_stats, index=tiles[tile_id_field])
        for stat in agg_stats:
            pattern = re.compile('tile_\d+_%s.tif' % stat)
            stat_match = [f.split('_')[1] for f in files if pattern.match(f)]
            tile_files[stat] = pd.Series(np.ones(len(stat_match)), index=stat_match)
        index_field = tiles.index.name
        tiles[index_field] = tiles.index
        tiles = tiles.set_index(tile_id_field, drop=False)[tile_files.isnull().any(axis=1)]
        tiles.set_index(index_field, inplace=True)#'''
    tiles['ul_x'] = [stem.get_ul_coord(xmin, xmax, x_res) 
                    for i, (xmin, xmax) in tiles[['xmin','xmax']].iterrows()]
    tiles['ul_y'] = [stem.get_ul_coord(ymin, ymax, y_res) 
                    for i, (ymin, ymax) in tiles[['ymin','ymax']].iterrows()]
    tiles['lr_x'] = [xmax if ulx == xmin else xmin for i, (ulx, xmin, xmax)
                    in tiles[['ul_x', 'xmin','xmin']].iterrows()]
    tiles['lr_y'] = [ymax if uly == ymin else ymin for i, (uly, ymin, ymax) 
                    in tiles[['ul_y', 'ymin','ymin']].iterrows()]
    
    support_nrows = int(support_size[0]/abs(y_res))
    support_ncols = int(support_size[1]/abs(x_res))
    t1 = time.time()
    #args = [(i + 1, n_tiles, t1, tile_info, mosaic_path, mosaic_tx, df_sets, df_var, (support_nrows, support_ncols), agg_stats, tile_path_template, prj, nodata, snap_coord) for i, (t_ind, tile_info) in enumerate(tiles[tiles['name'].isin(['1092', '3224'])].iterrows())]    
    args = [(i + 1, n_tiles, t1, tile_info, mosaic_path, mosaic_tx, df_sets, df_var, (support_nrows, support_ncols), agg_stats, tile_path_template, prj, nodata, snap_coord) for i, (t_ind, tile_info) in enumerate(tiles.iterrows())]
    
    if n_jobs > 1:
        print 'Predicting with %s jobs...\n' % n_jobs
        pool = Pool(n_jobs)
        limits = pool.map(stem.par_predict_tile, args, 1)
        pool.close()
        pool.join()
    else:
        print 'Predicting with 1 job ...\n'
        limits = []
        for arg in args:
            limits.append(stem.par_predict_tile(arg))#'''
    print '\n\nFinished predicting in %.1f hours. \n\nStitching tiles...' % ((time.time() - t1)/3600)

    limits = pd.concat(limits)
    t1 = time.time()
    mosaic_ul = mosaic_tx[0], mosaic_tx[3]
    driver = gdal.GetDriverByName('gtiff')
    for stat in agg_stats:
        dtype = mosaic.get_min_numpy_dtype(limits[stat])
        if stat == 'stdv':
            this_nodata = -9999
            ar = np.full((ysize, xsize), this_nodata, dtype=dtype) 
        else:
            this_nodata = nodata
            ar = np.full((ysize, xsize), this_nodata, dtype=dtype)
        
        for tile_id, tile_coords in tiles.iterrows():
            tile_file = os.path.join(tile_dir, 'tile_%s_%s.tif' % (tile_coords[tile_id_field], stat))
            ds = gdal.Open(tile_file)
            tile_tx = ds.GetGeoTransform()
            tile_ul = tile_tx[0], tile_tx[3]
            row_off, col_off = stem.calc_offset(mosaic_ul, tile_ul, mosaic_tx)
            # Make sure the tile doesn't exceed the size of ar
            tile_rows = min(ds.RasterYSize + row_off, ysize) - row_off
            tile_cols = min(ds.RasterXSize + col_off, xsize) - col_off
            ar_tile = ds.ReadAsArray(0, 0, tile_cols, tile_rows)
            try:
                ar[row_off : row_off + tile_rows, col_off : col_off + tile_cols] = ar_tile
            except Exception as e:
                import pdb; pdb.set_trace()
        
        out_path = os.path.join(out_dir, '%s_%s.tif' % (file_stamp, stat))
        #out_path = os.path.join('/home/server/pi/homes/shooper/delete_test', '%s_%s.tif' % (file_stamp, stat))
        gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(ar.dtype)
        mosaic.array_to_raster(ar, mosaic_tx, prj, driver, out_path, gdal_dtype, nodata=this_nodata)
    
    # Clean up the tiles
    shutil.rmtree(tile_dir)
    print 'Time for stitching: %.1f minutes\n' % ((time.time() - t1)/60)
    
    # Get feature importances and max importance per set
    t1 = time.time()
    print 'Getting importance values...'
    importance_cols = sorted([c for c in df_sets.columns if 'importance' in c])
    df_sets['max_importance'] = nodata
    if len(importance_cols) == 0:
        # Loop through and get importance
        importance_per_var = []
        for s, row in df_sets.iterrows():
            with open(row.dt_file, 'rb') as f: 
                dt_model = pickle.load(f)
            max_importance, this_importance = stem.get_max_importance(dt_model)
            df_sets.ix[s, 'max_importance'] = max_importance
            importance_per_var.append(this_importance)
        importance = np.array(importance_per_var).mean(axis=0)
    else:
        df_sets['max_importance'] = np.argmax(df_sets[importance_cols].values, axis=1)
        importance = df_sets[importance_cols].mean(axis=0).values
    pct_importance = importance / importance.sum()
    print '%.1f minutes\n' % ((time.time() - t1)/60)
    
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
            vote_path = os.path.join(out_dir, '%s_vote.tif' % file_stamp)
            ar_vote = gdal.Open(vote_path)
            print '\nComputing confusion matrix for vote...'
            vote_dir = os.path.join(model_dir, 'evaluation_vote')
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
            mean_path = os.path.join(out_dir, '%s_mean.tif' % file_stamp)
            ar_mean = gdal.Open(mean_path)
            print '\nGetting confusion matrix for mean...'
            mean_dir = os.path.join(model_dir, 'evaluation_mean')
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
    
    print '\nTotal prediction runtime: %.1f hours\n' % ((time.time() - t0)/3600)

if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))#'''