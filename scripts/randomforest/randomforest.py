import sys
import os
import time
import glob
import fnmatch
import gdal
import random
import pandas as pd
import numpy as np
import cPickle as pickle
from sklearn import ensemble
from sklearn import metrics
from matplotlib import pyplot as plt

import stem

def read_params(txt):
    '''
    Return a dictionary and a dataframe from parsed parameters in txt
    '''
    if not os.path.exists(txt):
        sys.exit('Param file does not exist:\n%s' % txt)

    d = {}

    # Read in the rest of the text file line by line
    try:
        with open(txt) as f:
            lines =  [line.split(';') for line in f]
            input_vars = [line.split(";") for line in f]
    except:
        print 'Problem reading parameter file:\n', txt
        return None

    # Set the dictionary key to whatever is left of the ";" and the val
    #   to whatever is to the right. Strip whitespace too.
    n_skip_lines = 0 #Keep track of the number of lines w/ a ";"
    #import pdb; pdb.set_trace()
    for var in lines:
        if len(var) == 2:
            d[var[0].strip()] = '"%s"' % var[1].strip().replace("\n", "")
            n_skip_lines += 1

    print '\nParameters read from:\n', txt, '\n'
    return d


def parse_bins(bin_str):
    ''' Integerize a string of min:max and return as a list of length 2'''
    bin_list = [b.split(':') for b in bin_str.split(',')]
    bins = [(int(mn), int(mx)) for mn, mx in bin_list]

    return bin


def random_values(ar, nodata, false_val, pct_train):
    """
    Return a tuple of arrays of unique values split randomly
    """
    vals = np.unique(ar[(ar != nodata) & (ar != false_val)])
    n_train = int(len(vals) * pct_train)
    train_vals = np.array(random.sample(vals, n_train))
    test_vals = vals[~np.in1d(vals, train_vals)]

    return train_vals, test_vals


def get_samples_from_zones(raster_path, col_name, data_band, n_samples, false_val, pct_train_sample, pct_train_zone, nodata=None, false_inflation=None):
    '''
    Return a dataframe of stratified randomly sampled pixels from raster_path
    '''
    if not os.path.exists(raster_path):
        sys.exit('Raster path specified does not exist: %s' % raster_path)
    print 'Reading the raster_path... %s\n%s\n' % (time.ctime(time.time()), raster_path)
    ds = gdal.Open(raster_path)
    tx = ds.GetGeoTransform()
    band = ds.GetRasterBand(data_band)
    ar = band.ReadAsArray()
    shape = ar.shape
    ds = None

    if nodata == None:
        nodata = band.GetNoDataValue()
        if nodata == None:
            sys.exit('Could not obtain NoData value from dataset and' + \
                     ' none specified in parameters file. Try re-running with' + \
                     'nodata specified.')

    samples_per = n_samples / 2

    print 'Making arrays of row and col indices... %s\n' % time.ctime(time.time())
    # Make 2D arrays of row and column indices
    ar_rows, ar_cols = np.indices(shape)

    # Get random samples for each bin
    train_rows = []
    train_cols = []
    test_rows = []
    test_cols = []
    train_zones, test_zones = random_values(ar, nodata, false_val, pct_train_zone)

    # Get the samples for the negative (condition is not present) values
    mask = ar == false_val
    false_rows = ar_rows[mask]
    false_cols = ar_cols[mask]
    try:
        if false_inflation:
            n_false_samples = int(samples_per * false_inflation)
            samples = random.sample(xrange(len(false_rows)), n_false_samples)
        else:
            n_false_samples = samples_per
            samples = random.sample(xrange(len(false_rows)), n_false_samples)
        #samples = random.sample(xrange(len(false_rows)), samples_per)
        these_rows = false_rows[samples]
        these_cols = false_cols[samples]
    # If there aren't enough pixels to generate samples for this bin
    except:
        print ('Not enough pixels equal to %s to generate %s' + \
               ' random samples. Returning all pixels for this bin.') \
              % (false_val, n_false_samples)
        these_rows = false_rows
        these_cols = false_cols
    # Split train and test for the false values
    split_ind = len(these_rows) * pct_train_sample
    train_rows.extend(these_rows[:split_ind])
    train_cols.extend(these_cols[:split_ind])
    test_rows.extend(these_rows[split_ind:])
    test_cols.extend(these_cols[split_ind:])

    # Get train and test samples for true values (condition is true)
    mask = np.in1d(ar, train_zones).reshape(shape)
    true_rows_train = ar_rows[mask]
    true_cols_train = ar_cols[mask]
    n_samples_train = int(samples_per * pct_train_sample)
    if true_rows_train.size <= n_samples_train:
        train_rows.extend(true_rows_train.ravel())
        train_cols.extend(true_cols_train.ravel())
    else:
        samples = random.sample(xrange(len(true_rows_train)), n_samples_train)
        train_rows.extend(true_rows_train[samples])
        train_cols.extend(true_cols_train[samples])
    #import pdb; pdb.set_trace()

    # Do the same for test samples of true values
    mask = np.in1d(ar, test_zones).reshape(shape)
    true_rows_test = ar_rows[mask]
    true_cols_test = ar_cols[mask]
    n_samples_test = int(samples_per * (1 - pct_train_sample))
    if true_rows_test.size <= n_samples_test:
        test_rows.extend(true_rows_test.ravel())
        test_cols.extend(true_cols_test.ravel())
    else:
        samples = random.sample(xrange(len(true_rows_test)), n_samples_test)
        test_rows.extend(true_rows_test[samples])
        test_cols.extend(true_cols_test[samples])
    #import pdb; pdb.set_trace()

    # Calculate x and y for later extractions. Then get array values.
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    train_x = [int(ul_x + c * x_res) for c in train_cols]
    train_y = [int(ul_y + r * y_res) for r in train_rows]
    ar[np.in1d(ar, np.concatenate((train_zones, test_zones))).reshape(shape)] = 1 #Set any zone values to 1
    train_vals = ar[train_rows, train_cols]
    df_train = pd.DataFrame(zip(train_x, train_y, train_rows, train_cols, train_vals),
                            columns=['x', 'y', 'row', 'col', col_name])
    test_x = [int(ul_x + c * x_res) for c in test_cols]
    test_y = [int(ul_y + r * y_res) for r in test_rows]
    test_vals = ar[test_rows, test_cols]
    df_test = pd.DataFrame(zip(test_x, test_y, test_rows, test_cols, test_vals),
                           columns=['x', 'y', 'row', 'col', col_name])

    return df_train, df_test, x_res



def get_stratified_sample(raster_path, col_name, data_band, n_samples, bins, pct_train=None, nodata=None):
    '''
    Return a dataframe of stratified randomly sampled pixels from raster_path
    '''
    if not os.path.exists(raster_path):
        sys.exit('Raster path specified does not exist: %s' % raster_path)
    print 'Reading the raster_path... %s\n%s\n' % (time.ctime(time.time()), raster_path)
    ds = gdal.Open(raster_path)
    tx = ds.GetGeoTransform()
    band = ds.GetRasterBand(data_band)
    ar = band.ReadAsArray()
    ds = None

    if nodata == None:
        nodata = band.GetNoDataValue()
        if nodata == None:
            sys.exit('Could not obtain NoData value from dataset and' + \
                     ' none specified in parameters file. Try re-running with' + \
                     'nodata specified.')

    samples_per = n_samples / len(bins)

    print 'Making arrays of row and col indices... %s\n' % time.ctime(time.time())
    # Make 2D arrays of row and column indices
    ar_rows, ar_cols = np.indices(ar.shape)

    # Get random samples for each bin
    train_rows = []
    train_cols = []
    test_rows = []
    test_cols = []
    #nodata_mask = ar != nodata
    for b in bins:
        this_min, this_max = b
        print 'Getting random samples between %s and %s...' % (this_min, this_max)
        print time.ctime(time.time()), '\n'
        mask = (ar > this_min) & (ar <= this_max) & (ar != nodata)
        these_rows = ar_rows[mask]
        these_cols = ar_cols[mask]

        try:
            samples = random.sample(xrange(len(these_rows)), samples_per)
            tr_rows = these_rows[samples]
            tr_cols = these_cols[samples]
        # If there aren't enough pixels to generate samples for this bin
        except:
            print ('Not enough pixels between %s and %s to generate %s' + \
                   ' random samples. Returning all pixels for this bin.') \
                  % (this_min, this_max, samples_per)
            tr_rows = these_rows
            tr_cols = these_cols

        # If pct_train is specified, split the sample indices into train/test sets
        te_rows = []
        te_cols = []
        if pct_train:
            split_ind = len(tr_rows) * pct_train
            te_rows = tr_rows[split_ind:]
            te_cols = tr_cols[split_ind:]
            tr_rows = tr_rows[:split_ind]
            tr_cols = tr_cols[:split_ind]

        train_rows.extend(tr_rows)
        train_cols.extend(tr_cols)
        test_rows.extend(te_rows)
        test_cols.extend(te_cols)

    # Calculate x and y for later extractions. Then get array values.
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx

    train_x = [int(ul_x + c * x_res) for c in train_cols]
    train_y = [int(ul_y + r * y_res) for r in train_rows]
    train_vals = ar[train_rows, train_cols]
    df_train = pd.DataFrame(zip(train_x, train_y, train_rows, train_cols, train_vals),
                            columns=['x', 'y', 'row', 'col', col_name])
    df_test = None
    if pct_train:
        test_x = [int(ul_x + c * x_res) for c in test_cols]
        test_y = [int(ul_y + r * y_res) for r in test_rows]
        test_vals = ar[test_rows, test_cols]
        df_test = pd.DataFrame(zip(test_x, test_y, test_rows, test_cols, test_vals),
                               columns=['x', 'y', 'row', 'col', col_name])

    return df_train, df_test, x_res


def extract_rowcol(df, predictor_name, raster_path):
    ''' Sample the predictor values at each row/col pair in df '''

    ds = gdal.Open(raster_path)
    tx = ds.GetGeoTransform()
    ar = ds.ReadAsArray()
    ds = None

    columns = df.columns
    if 'row' not in columns and 'col' not in columns:
        sys.exit('KeyError: Either "row" or "col" is not a column in the passed dataframe')
    df[predictor_name] = ar[df.row, df.col]


def sample_predictors(df, df_var, nodata):
    ''' Add samples in place for each predictor in df_var to df from [df.row, df.col]'''

    t0 = time.time()
    print 'Sampling predictors...'
    for var, row in df_var.iterrows():
        print 'Getting samples for ', var
        this_path = row.file
        if not os.path.exists(this_path):
            sys.exit('\nERROR: Raster path specified does not exist: %s' % this_path)
        ds = gdal.Open(this_path)
        this_nodata = row.nodata
        band = ds.GetRasterBand(row.databand)
        ar = band.ReadAsArray()
        ds = None
        df[var] = ar[df.row, df.col]
        df.ix[df[var] == this_nodata, var] = nodata
    print 'Time for sampling predictors: %.1f seconds\n' % (time.time() - t0)

    return df[df[var] != nodata]


def train_rf_classifier(x, y, ntrees=50, njobs=12, max_depth=None, max_features='auto'):
    ''' Return a trained Random Forest classifier'''
    rf = ensemble.RandomForestClassifier(n_estimators=ntrees,
                                         n_jobs=njobs,
                                         max_depth=max_depth,
                                         max_features=max_features,
                                         oob_score=True)
    rf.fit(x, y)

    return rf


def train_rf_regressor(x, y, ntrees=50, njobs=12, max_depth=None, max_features='sqrt'):
    ''' Return a trained Random Forest regressor'''
    rf = ensemble.RandomForestRegressor(n_estimators=ntrees,
                                        n_jobs=njobs,
                                        max_depth=max_depth,
                                        max_features=max_features,
                                        oob_score=True)
    rf.fit(x, y)

    return rf

def save_rfmodel(rf, filename):
    ''' Write a RandomForest model to disk'''

    with open(filename, 'w+') as f:
        pickle.dump(rf, f, protocol=-1)

    return filename


def find_file(basepath, search_str, tsa_str=None, path_filter=None):
    '''
    Return the full path within the directory tree /basepath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    '''if not os.path.exists(basepath):
        print 'basepath does not exist: \n%s' % basepath
        return None'''
     
    if tsa_str: 
        bp = os.path.join(basepath, tsa_str)

        # Search the bp directory tree. If search_str is in a file, get the full path.
        paths = []
        for root, dirs, files in os.walk(bp, followlinks=True):
            these_paths = [os.path.join(root, f) for f in files]
            these_paths = fnmatch.filter(these_paths, search_str)
            paths.extend(these_paths)
    else:
        paths = glob.glob(os.path.join(basepath, search_str))
    
    # If path filter is specified, remove any paths that contain it
    if path_filter:
        [paths.remove(p) for p in paths if fnmatch.fnmatch(p, path_filter)]
        #paths = [p for p in paths if fnmatch.fnmatch(p, path_filter)]
    
    '''if len(paths) > 1:
        print 'Multiple files found for tsa: ' + tsa_str
        for p in paths:
            print p
        print 'Selecting the first one found...\n'# '''
    
    #import pdb; pdb.set_trace()
    if len(paths) == 0:
        #pdb.set_trace()
        sys.exit(('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str))
    
    return paths[0]



def get_predictors(df_var, nodata, ydims=None, constant_vars=None):
    '''
    Return an array of flattened predictor arrays where each predictor is a
    separate column
    '''
    t0 = time.time()
    predictors = []
    for var, row in df_var.iterrows():
        data_band, search_str, basepath, by_tsa, path_filter = row
        #data_band, search_str, basepath, by_tsa, nodata, path_filter = row
        
        if constant_vars: search_str = search_str.format(constant_vars['YEAR'])        
        
        this_path = find_file(basepath, search_str)
        if not os.path.exists(this_path):
            sys.exit('\nERROR: Raster path specified does not exist: %s' % this_path)
        ds = gdal.Open(this_path)
        try:
            xsize = ds.RasterXSize
        except:
            import pdb; pdb.set_trace()
        if ydims:
            upper_ydim, lower_ydim = ydims
            this_ysize = lower_ydim - upper_ydim
            ar_var = ds.ReadAsArray(0, upper_ydim, xsize, this_ysize)
        else:
            ar_var = ds.ReadAsArray()
        ds = None
        if 'nodata' in row:
            this_nodata = row.nodata
            ar_var[ar_var == this_nodata] = nodata # Change var nodata to nodata val for output
        predictors.append(ar_var.ravel())
        #print 'Shape of %s: %s' % (var, ar_var.ravel().shape)
    
    if constant_vars:
        size = predictors[0].size
        for const in sorted(constant_vars.keys()):
            val = constant_vars[const]
            predictors.append(np.full(size, val, dtype=np.int16))
            
    # Make an array where each row is a different predictor
    ar = np.vstack(predictors).T
    del predictors
    nodata_mask = np.all(~(ar == nodata), axis=1) #Maybe should be np.any?
    ar = ar[nodata_mask]

    print 'Finished getting arrays: %.1f minutes' % ((time.time() - t0)/60)
    return ar, nodata_mask

def par_predict_from_dt(args):
    dt, predictors = args
    return dt.predict(predictors)

def array_to_raster(array, tx, prj, driver, out_path, dtype, nodata=None):
    ''' Save a numpy array as a new raster '''
    print 'Saving raster...'
    rows,cols = array.shape
    #save new raster
    out_ds = driver.Create(out_path, cols, rows, 1, dtype)
    if out_ds is None:
        print sys.exit('\nCould not create ' + out_path)
    #write the data
    band = out_ds.GetRasterBand(1)
    band.WriteArray(array)
    #flush data to disk
    band.FlushCache()
    if nodata != None: band.SetNoDataValue(nodata)

    #georeference the image and set the projection
    out_ds.SetGeoTransform(tx)
    out_ds.SetProjection(prj)
    print 'Raster written to:\n', out_path


def calc_rmse(ar, test_samples, val_col, nodata):
    '''
    Return the Root Mean Squared Error for predicted values from ar
    and observed values from test_samples.
    '''
    #test_samples = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    p_vals = ar[test_samples.row, test_samples.col]
    mask = (p_vals != nodata)
    p_vals = p_vals[mask]
    t_vals = test_samples.ix[mask, val_col].values#'''

    rmse  = round(np.sqrt(((t_vals - p_vals) ** 2)).mean(), 3)
    t_mask = t_vals == 1 # Where observed value is true (is a landslide)
    rmse_true  = round(np.sqrt(((t_vals[t_mask] - p_vals[t_mask]) ** 2)).mean(), 3)
    rmse_false = round(np.sqrt(((t_vals[~t_mask] - p_vals[~t_mask]) ** 2)).mean(), 3)

    return rmse, rmse_true, rmse_false


def calc_auc(ar, test_samples, val_col, nodata, out_dir=None):
    '''
    Return the AUC score for predicted values from ar and observed values
    from test_samples. Optionally, plot the ROC curve if out_dir is
    specified.
    '''
    p_vals = ar[test_samples.row, test_samples.col]
    mask = (p_vals != nodata)
    p_vals = p_vals[mask]
    t_vals = test_samples.ix[mask, val_col]

    auc = round(metrics.roc_auc_score(t_vals, p_vals), 3)
    if out_dir:
        fpr, tpr, _ = metrics.roc_curve(t_vals, p_vals)
        plt.plot(fpr, tpr, '-')
        plt.axis([0, 1, 0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        out_png = os.path.join(out_dir, 'roc_curve.png')
        plt.savefig(out_png)
        print 'AUC plot written to: ', out_png
        plt.clf()

    return auc


def calc_brier(ar, sample_txt, val_col):

    test_samples = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    p_vals = ar[test_samples.row, test_samples.col]
    mask = p_vals != nodata
    p_vals = p_vals[mask]
    t_vals = test_samples.ix[mask, val_col]

    brier_score = ((t_vals - p_vals) ** 2).mean()

    return brier_score
















