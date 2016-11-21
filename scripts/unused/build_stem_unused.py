# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:29:36 2016

@author: shooper
"""

    c = 1
    
    
    '''for ind, row in df_sets.iterrows():
        print 'Getting arrays for set %s of %s' % (c, len(support_sets))
        ar_coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
        tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, ar_coords, mosaic_tx, xsize, ysize)
        # Get the ids of TSAs this kernel covers
        tsa_ids = np.unique(tsa_ar)
        tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=0]
        print tsa_strs
        
        pred_arrays = []
        for var in predict_cols:
            files = [f for f in np.unique(df_train['file_' + var]) for tsa in tsa_strs if tsa in f]
            print files
            
            pred_arrays.append(mosaic.get_mosaic(mosaic_tx, tsa_strs, tsa_ar, tsa_off, ar_coords, tsa_txt,  data_band, files=files).ravel())
        df_sets.ix[ind, predict_cols] = pred_arrays
        c+=1'''
    
    print '\n', time.time() - t0
    # Make prediction
    #x_test = df_ts[predict_cols]
    #predictions = dt.predict(x_test)
    '''df_pred = pd.DataFrame(columns = predict_cols)
    static_cols = ['YEAR', 'DAY', 'TIME', 'EFFORT_HRS']
    val_cols = [predict_cols.remove(i) for i in static_cols]
    files = glob.glob(os.path.join(search_dir, '*') + search_str)
    for f in files:
        bn = os.path.basename(f).replace('.tif', '')
        ar, shape = raster_to_1darray(f)
        df_pred[bn] = ar
    
    df_pred[static_cols] = np.array([2012, 150, 8, 2] * len(df_pred)).reshape(len(df_pred), 4)
    cols = list(df.columns)
    [cols.remove(i) for i in ['x', 'y', target_col]]
    df_pred = df_pred[cols]
    pred_array = dt.predict(df_pred).reshape(shape)'''
    
def get_predict_arrays(df_train, predict_cols, mosaic_tx, tsa_strs, tsa_ar, tsa_off, ar_coords,  data_band):
    '''
    Return a dataframe where each columns is composed of a 1D array from
    mosaic_by_tsa for each predictor variable in predict_cols
    '''
    df_predict = pd.DataFrame(columns=predict_cols)
    n_vars = len(predict_cols)
    for i, var in enumerate(predict_cols):
        print 'Getting array for %s of %s variables' % (i+1, n_vars)
        files = [f for f in np.unique(df_train['file_' + var])
        for tsa in tsa_strs if tsa in f]
        
        df_predict[var] = mosaic.get_mosaic(mosaic_tx, tsa_strs, tsa_ar,\
        tsa_off, ar_coords,  data_band, files=files).ravel()
        
    return df_predict


def predict_set(set_id, df_train, predict_cols, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, data_band, tree):
    '''
    Return a predicted array for set, set_ind
    '''
    #ar_coords = df_sets.ix[set_id, ['ul_x', 'ul_y', 'lr_x', 'lr_y']]
    
    # Get an array of tsa_ids within the bounds of ar_coords
    tsa_ar, tsa_off = mosaic.extract_kernel(mosaic_ds, ar_coords,
                                            mosaic_tx, xsize, ysize)
    # Get the ids of TSAs this kernel covers
    tsa_ids = np.unique(tsa_ar)
    tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=0]
    print tsa_strs
    
    df_predict = get_predict_arrays(df_train, predict_cols, mosaic_tx, tsa_strs, tsa_ar,
                                    tsa_off, ar_coords, data_band)
    #tree = df_sets.ix[set_id, 'tree']
    prediction = tree.predict(df_predict).reshape(tsa_ar.shape)
    
    return set_id, prediction





#    set_id = support_sets[0]
#    ar_coords = df_sets.ix[set_id, ['ul_x', 'ul_y', 'lr_x', 'lr_y']]
#    tree = df_sets.ix[set_id, 'tree']

#set_id, ar = predict_set(set_id, df_train, predict_cols, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, data_band, tree)

#file_dict = {var: [f for f in np.unique(df_train.ix[df_train.set_id==set_id, 'file_' + var])] for var in predict_cols}

#   First, get anything with 'file' in the name (i for i in df loops through columns)
file_cols = [col for col in df_train if 'file' in col] 

# Then get all columns if the column name is in a string in file_cols
#    but the columns name itself is not in file_cols
predict_cols = list(np.unique([col_p for col_p in df_train for col_f in file_cols
if col_p in col_f and col_p not in file_cols])):
    
    
''' from get_overlapping_sets'''
if tx[1] > 0:
    min_x = 'ul_x'
    max_x = 'lr_x'
else:
    min_x = 'lr_x'
    max_x = 'ul_x'
if tx[5] > 0:
    min_y = 'ul_y'
    max_y = 'lr_y'
else:
    min_y = 'lr_y'
    max_y = 'ul_y'

df_overlap = df_sets[
(df_sets[min_x] > row[min_x]) & (df_sets[min_x] < row[max_x]) |
(df_sets[max_x] > row[min_x]) & (df_sets[max_x] < row[max_x]) |
(df_sets[min_y] > row[min_y]) & (df_sets[min_y] < row[max_y]) |
(df_sets[max_y] > row[min_y]) & (df_sets[max_y] < row[max_y]) ]

def parallelize_predict_set(args):
    '''Helper function to parallelize predict_set()'''
    try:
        return predict_set(*args)
    except:
        print 'Problem with set_id ', args[0]
        traceback.print_exception(*sys.exc_info())
        print '\n'

def get_predict_arrays(df_var, mosaic_tx, tsa_strs, tsa_ar, tsa_off, ar_coords):
    '''
    Return a dataframe where each columns is composed of a 1D array from
    mosaic_by_tsa for each predictor variable in predict_cols
    '''
    df_predict = pd.DataFrame(columns=df_var.index)
    n_vars = len(df_var)
    for i, var in enumerate(df_var.index):
        print 'Getting array for %s of %s variables' % (i+1, n_vars)
        data_band, search_str, basepath, path_filter = df_var.ix[var]
        
        files = [find_file(basepath, tsa, search_str, path_filter) for tsa in tsa_strs]
        
        df_predict[var] = mosaic.get_mosaic(mosaic_tx, tsa_strs, tsa_ar,\
         ar_coords,  data_band, files=files).ravel()
        
    return df_predict

def par_predict_set(args):
    ''' Helper function for prediction for a decision tree'''
    inds, df, dt = args
    prediction = dt.predict(df)
    return inds, prediction.astype(np.int32)

def pieces_of_predict_set():   
    # Split the predictor df into chunks, and predict for each chunk
    #   First get indices for each chunk
    '''n_pixels = len(df_predict)
    chunk_size = n_pixels/44 #44 is number of processors to use
    inds_left = np.array(range(0, n_pixels, chunk_size)[:-1])
    inds_right = inds_left + chunk_size
    inds_right[-1] = n_pixels - 1
    inds = zip(inds_left, inds_right)  
    
    # Predict in parallel
    t0 = time.time()
    par_args = [[i, df_predict.ix[i[0] : i[1] - 1], dt] for i in inds]
    p = Pool(44)
    prediction_iter = p.imap(par_predict_set, par_args)
    prediction_pieces = [it for it in prediction_iter]
    p.close()
    p.join()
    del df_predict, par_args# Release predictor array resources
    t1 = time.time()'''

    '''# Stitch the pieces back together
    ar_prediction = np.empty(n_pixels, dtype=np.int32)
    for piece in prediction_pieces:
        inds, ar = piece
        ar_prediction[inds[0]:inds[1]] = ar
        del ar, piece
    del prediction_pieces
    ar_prediction = ar_prediction.reshape(array_shape)
    print 'Finished stitching...', time.time() - t1
    print 'Data type: ', ar_prediction.dtype'''
    
''' from getting arrays'''
        #predictors.append((var, ar_var.ravel()))
        #args = [var, mosaic_tx, tsa_strs, tsa_ar, tsa_off, ar_coords, data_band, files]
        #pred_args.append(args)
        #predictors.append(get_predict_array(args))
    
    '''t0 = time.time()
    p = Pool(n_workers)
    predictors_iter = p.imap(get_predict_array, pred_args)
    predictors = [it for it in predictors_iter]
    p.close()
    p.join()'''
    #df_predict = pd.DataFrame({v: ar for v, ar in predictors})


def par_fill_tile_band(args):
    '''Helper function to be able to run fill_tile_band() in parallel'''
    return fill_tile_band(*args)
    


'''for set_id, ar in predictions.iteritems():
    
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    ul_x, ul_y = df_sets.ix[set_id, ['ul_x', 'ul_y']]
    ul_x += (ul_x - m_ulx) % x_res
    ul_y += (ul_y - m_uly) % y_res
    tx = ul_x, x_res, x_rot, ul_y, y_rot, y_res
    out_path = out_dir + '/dt_test_%s.bsq' % set_id
    mosaic.array_to_raster(ar, tx, prj, driver, out_path, GDT_Int32, nodata=nodata)'''

'''# For testing aggregation without generating predictions
predictions = {}
for i in range(n_sets):
    this_ar = np.full((13333,10000), nodata, dtype=np.int32)
    r = random.sample(xrange(13333), 2)
    c = random.sample(xrange(10000), 2)
    ar_rand = np.random.randint(0, 100, (max(r) - min(r), max(c) - min(c)))
    this_ar[min(r):max(r), min(c):max(c)] = ar_rand
    predictions[list(df_these_sets.index)[i]] = this_ar#'''

def predict_set_from_disk(df_sets, set_id, params):
    
    inputs, df_var = read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    df_var = df_var.reindex(df_var.index.sort_values())
    this_set = df_sets.ix[set_id]
    with open(this_set.dt_file, 'rb') as f: 
        dt_model = pickle.load(f)
    
    mosaic_ds = gdal.Open(mosaic_path)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    xsize = mosaic_ds.RasterXSize
    ysize = mosaic_ds.RasterYSize
    prj = mosaic_ds.GetProjection()
    driver = mosaic_ds.GetDriver()
    
    ar_coords = this_set[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
    mosaic_dir = '/vol/v2/stem/canopy/canopy_20160212_2016/var_mosaics'
    saving_stuff = set_id, mosaic_dir, prj, driver
    ar_predict = predict_set(set_id, df_var, mosaic_ds, ar_coords, mosaic_tx, xsize, ysize, dt_model, saving_stuff)
    return ar_predict
    
    '''out_dir = '/vol/v2/stem/scripts/testing'
    out_path = os.path.join(out_dir, 'predict_rerun_%s.bsq' % set_id)
    
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    tx = this_set.ul_x, x_res, x_rot, this_set.ul_y, y_rot, y_res
    mosaic.array_to_raster(ar_predict, tx, prj, driver, out_path, GDT_Int32)'''

def main(params):
    
    '''### copy params to out_dir #### '''
    
    #read_params(params)
    inputs, df_var = read_params(params)

    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    try:
        num_vars = vars_to_numbers(cell_size, support_size, sets_per_cell,
                                   min_obs, pct_train, n_tiles, nodata)
        cell_size, support_size, sets_per_cell, min_obs, pct_train, n_tiles, nodata = num_vars
        str_check = sample_txt, target_col, mosaic_path, tsa_txt, out_dir
    except NameError as e:
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)
        return None
    
    now = datetime.now()
    date_str = str(now.date()).replace('-','')
    time_str = str(now.time()).replace(':','')[:4]
    stamp = '{0}_{1}_{2}'.format(target_col, date_str, time_str)
    out_dir = os.path.join(out_dir, stamp)
    os.makedirs(out_dir) # With a timestamp in dir, no need to check if it exists
    shutil.copy2(params, out_dir) #Copy the params for reference
    
    # Get samples and support set bounds
    if 'gsrd_shp' not in locals(): gsrd_shp = None
    out_txt = os.path.join(out_dir, stamp + '.txt')
    dfs = gsrd.get_gsrd(mosaic_path, cell_size, support_size, sets_per_cell,
                        sample_txt, min_obs, pct_train, target_col, out_txt,
                        gsrd_shp)
    df_train, df_test, df_sets = dfs
    support_sets = df_train.set_id.unique()

    # Check that df_train has exactly the same columns as variables specified in df_vars
    unmatched_vars = [v for v in df_var.index if v not in [c for c  in df_train]]
    
    if len(unmatched_vars) != 0:
        unmatched_str = '\n'.join(unmatched_vars)
        msg = 'Columns not in sample_txt but specified in params:\n' + unmatched_str
        raise NameError(msg)
    
    predict_cols = sorted(np.unique([c for c in df_train.columns for v in df_var.index if v in c]))
    df_var = df_var.reindex(df_var.index.sort_values())# Make sure predict_cols and df_var are in the same order

    # Train a tree for each support set
    x_train = df_train.reindex(columns=predict_cols + ['set_id'])
    y_train = df_train[[target_col, 'set_id']]    
    df_sets['dt_model'] = [fit_tree(x_train.ix[x_train.set_id==s, predict_cols],\
    y_train.ix[y_train.set_id==s, target_col]) for s in support_sets]
    
    # Write df_sets and each decison tree to disk
    write_model(out_dir, df_sets)
    
    
    ''' Split script here to be able to predict for other years'''
    mosaic_ds = gdal.Open(mosaic_path, GA_ReadOnly)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    xsize = mosaic_ds.RasterXSize
    ysize = mosaic_ds.RasterYSize
    prj = mosaic_ds.GetProjection()
    driver = mosaic_ds.GetDriver()
    
    t0 = time.time()
    
    predict_dir = os.path.join(out_dir, 'predctions')
    os.mkdir(predict_dir)
    # Loop through each set and generate predictions
    m_ulx, x_res, x_rot, m_uly, y_rot, y_res = mosaic_tx
    c = 1
    total_sets = len(support_sets)
    predictions = {}
    for set_id, row in df_sets.iterrows():
        print 'Predicting for set %s of %s' % (c, total_sets)
        ar_coords = row[['ul_x', 'ul_y', 'lr_x', 'lr_y']]
        ar_predict = predict_set(set_id, df_var, mosaic_ds, ar_coords, 
                                 mosaic_tx, xsize, ysize, row.dt_model, nodata)
        #predictions[set_id] = ar_predict
        
        tx = ar_coords['ul_x'], x_res, x_rot, ar_coords['ul_y'], y_rot, y_res
        out_path = predict_dir + '/prediction_%s.bsq' % set_id
        mosaic.array_to_raster(ar_predict, tx, prj, driver, out_path, GDT_Int32, nodata=nodata)
        c += 1
    mosaic_ds = None                  
    print '\nTotal time for predictions: %.1f minutes' % ((time.time() - t0)/60)#'''
    
    #Aggregate predictions by tile and stitch them back together
    aggr.aggregate_predictions(ysize, xsize, nodata, n_tiles, mosaic_tx, support_size, predict_dir, df_sets, out_dir, stamp, prj, driver)
    

if __name__ == '__main__':
     params = sys.argv[1]
     sys.exit(main(params))