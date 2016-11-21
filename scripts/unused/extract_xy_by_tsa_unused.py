# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:22:08 2016

@author: shooper
"""

""" Attempts at extracting that take way too long"""
# Get data values at each row,col
# Group by raster cell
print 'extracting...'
print time.ctime(time.time())
#gb = df_xy.groupby('rowcol') 
gb = np.unique(df_xy.loc[df_xy['file'] == f, 'rowcol'])
# Make a dictionary of {rowcol: data val}
val_dict = dict([(k, extract_by_rowcol(k, ar, nodata)) 
for k in gb])#.groups.keys()])
for k in val_dict:
    df_xy.loc[df_xy['rowcol'] == k, 'val'] = val_dict[k]
# Make a list of data values of the same
#l = [vals.extend([v] * len(df[df['rowcol'] == k])) for k, v in val_dict.items()]
#df_xy.loc[df_xy['file'] == f, 'val'] = df_xy['rowcol'].apply(
#lambda rc: extract_by_rowcol(rc, ar, nodata))
#df_xy.loc[df_xy['file'] == f, 'val'] = df_xy[['x', 'y']].apply(
#lambda xy: extract_at_xy(xy, ar, tx, nodata, val_only=True), axis=1)
#df_xy.loc[df_xy['file'] == f, 'rowcol'] = df_xy[['row', 'col']]\
#.apply(lambda rc: (rc[0] - row_off, rc[1] - col_off), axis=1)
    
def extract_by_rowcol(rowcol, array, nodata):
    
    row, col = rowcol

    try:
        val = array[row, col]
    except IndexError:
        val = nodata

    return val

#df_temp['val'] = ar[df_temp['data_row'], df_temp['data_col']]
#df_temp[val_cols] = np.array([[ar[r + dr, c + dc] for dr, dc in kernel_dirs] for r, c in df_temp[['data_row', 'data_col']]])#.reshape(len(val_cols, len))

        '''print 'Extracting for file ', f
        print time.ctime(time.time()), '\n'
        ds = gdal.Open(f, GA_ReadOnly)
        tx = ds.GetGeoTransform()
        ar = ds.GetRasterBand(data_band).ReadAsArray()
        
        #Calc row, col of the dataset for each point
        row_off = int((tx[3] - mosaic_tx[3])/tx[5])
        col_off = int((tx[0] - mosaic_tx[0])/tx[1])
        data_rows = [row - row_off + d for row in df_temp['row'] for d in row_dirs]
        data_cols = [col - col_off + d for col in df_temp['col'] for d in col_dirs]
        
        # Get the data values
        vals = ar[data_rows, data_cols].reshape(len(df_temp), len(val_cols))
        df_vals = pd.DataFrame(vals, columns=val_cols, index=df_temp.index)
        df_temp[val_cols] = df_vals
        dfs.append(df_temp)
        
        ds = None
        # Reassign the temp df back to the main df
        #df_xy[df_xy['file'] == f] = df_temp'''

    '''df_xy[['row', 'col', 'tsa_id']] = df_xy[['x', 'y']].apply(
    lambda xy: extract_at_xy(xy, mosaic_ar, mosaic_tx, nodata), axis=1)
    df_xy = df_xy[df_xy['tsa_id'] != nodata]#'''

def calc_offset(xy, data_tx):
    '''
    Return the row and col of xy in a dataset given a geotransform
    '''
    x, y = xy
    ul_x, ul_y = data_tx[0], data_tx[3]
    
    row = int((y - ul_y)/data_tx[5])
    col = int((x - ul_x)/data_tx[1])
    
    return pd.Series((row, col))

def extract_at_xy(xy, array, data_tx, nodata, val_only=False):
    '''
    Return the pixel value at x,y in the array
    '''
    x, y = xy
    ul_x, ul_y = data_tx[0], data_tx[3]
    
    row = int((y - ul_y)/data_tx[5])
    col = int((x - ul_x)/data_tx[1])
    
    try:
        val = array[row, col]
    except IndexError:
        val = nodata
    
    if val_only:
        return val
    else:
        return pd.Series((row, col, val))
        
def get_random_pixels(min_val, max_val, extent_array, tx, n_xbins, n_ybins, n_samples):
    '''
    Return stratified random samples of x,y coords within the extent of
    the dataset
    '''
    # Get specs from the dataset
    x_size = ds.RasterXSize
    y_size = ds.RasterYSize
    ul_x, x_res, x_rot, ul_y, y_rot, y_res = tx
    lr_x = ul_x + x_size * x_res
    lr_y = ul_y + y_size * y_res
    
    # Get min and max x and y values
    min_x = min(ul_x, lr_x)
    min_y = min(ul_y, lr_y)
    max_x = max(ul_x, lr_x)
    max_y = max(ul_y, lr_y)
    
    # Calculate the number of samples per bin and width of bins
    x_range = max_x - min_x
    y_range = max_y - min_y
    samples_per = n_samples/(n_xbins * n_ybins)
    x_step = x_range/n_xbins #increase per bin
    y_step = y_range/n_ybins
    x_mins = [min_x + i * x_step for i in xrange(n_xbins)]
    y_mins = [min_y + i * y_step for i in xrange(n_ybins)]
    
    samples = []
    i = 0
    while i < n_bins:
        this_max_x = min_x + x_step
        this_max_y = 
        x = random.sample(xrange(min_x, max_x + 1), n_samples)
        y = random.sample(xrange(min_y, max_y + 1), n_samples)

# If the years are represented in different bands, get the band
#  from the metada
'''if data_band < 0:
    band_dict = ds.GetMetadata_Dict()
    inverse = {int(v): int(k.split('_')[1]) for k,v in band_dict.iteritems()}
    data_band = inverse[year]'''
    
'''def main(params, out_dir=None, xy_txt=None):          
    
    # Read params. Make variables from each line of the 1-line variables
    inputs, df_vars = read_params(params)
    for var in inputs:
        exec ("{0} = str({1})").format(var, inputs[var])
    
    if not os.path.exists(out_dir):
        print 'Output directory does not exist:\n', out_dir
        return None
    
    if 'years' in locals(): 
        years =  [int(yr) for yr in years.split(',')]
    else:
        try:
            year_start = int(year_start)
            year_end = int(year_end)
            years = range(year_start, year_end + 1)
        except NameError:
            print 'No list of years or year_start/year_end specified in' +\
            ' param file:\n%s\n. Re-run script with either of these' +\
            ' parameters set.' % params
            return None
    
    # Get the TSA mosaic as an array
    print 'Reading mosaic dataset...%s\n' % time.ctime(time.time())
    mosaic_ds = gdal.Open(tsa_mosaic, GA_ReadOnly)
    mosaic_tx = mosaic_ds.GetGeoTransform()
    mosaic_ar = mosaic_ds.ReadAsArray()
    
    # Read in the TSA txt and XY txt as dataframes
    print 'Reading in text files... %s\n' % time.ctime(time.time())
    df_tsa = pd.read_csv(tsa_txt, sep='\t', dtype={'tsa_str':object})
    df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')
    #if xy_txt: 
    #    df_xy = pd.read_csv(xy_txt, sep='\t', index_col='obs_id')
    
    # Only need to do this once per xy sample set
    # Get the TSA at each xy location and the row and col 
    print 'Extracting tsa array values... %s\n' % time.ctime(time.time())
    extract_at_xy(df_xy, mosaic_ar, mosaic_tx, 'tsa_id')# Gets val inplace
    
    # For each year, do extractions
    tsa_ids = np.unique(df_tsa.tsa_id) # Assumes there are coords in all TSAs
    c = 1 #For keeping count of files processed
    n_years = len(years)
    n_vars = len(df_vars)
    n_files = len(df_tsa) * n_years * n_vars
    xy_txt_bn = os.path.basename(xy_txt)
    out_dir = os.path.join(out_dir, xy_txt_bn.split('.')[0])
    if not os.path.exists(out_dir): os.mkdir(out_dir)
    
    for var_name, var_row in df_vars.iterrows(): #var_name is index col
        # Get variables from row
        search_str  = var_row.search_str
        basepath    = var_row.basepath
        path_filter = var_row.path_filter
        data_type   = var_row.data_type
        data_band   = var_row.data_band
        
        df_var = pd.DataFrame()
        for year in years:
            if var_row.data_band < 0: 
                data_band = year - 1984 + 1
            else: 
                data_band = var_row.data_band
            
            dfs = [] # For storing kernel and stats
            print 'Finding data files and extracting for %s... %s\n' % (year, datetime.now())
            var_col = var_name + str(year)
            file_col = 'file_' + var_col
            # Store the filepath for each TSA
            df_tsa[file_col] = [find_file(basepath, tsa, search_str.format(year), path_filter)
            for tsa in df_tsa['tsa_str']] 
            
            # Handle any rows for for which the file is null
            if df_tsa[file_col].isnull().any():
                df_null = df_tsa[df_tsa[file_col].isnull()]
                print 'TSAs excluded from extractions for %s from %s:' % (var_name, year)
                for ind, row in df_null.iterrows(): print row['tsa_str']
                print ''
                n_null = len(df_null)
                # Make the file name a unique integer so that it can be
                #   distinguished from real files and from other null files
                df_tsa.loc[df_null.index, file_col] = range(n_null)
                
            # Get the file string for each xy for each year
            df_xy[file_col] = ''# Creates the column but keeps it empty
            for tsa in tsa_ids:
                df_xy.loc[df_xy['tsa_id'] == tsa, file_col] = df_tsa.loc[
                df_tsa['tsa_id'] == tsa, file_col].values[0]
        
            # For each file, get the dataset as an array and extract all values at each row col
            val_cols = ['%s_%s' % (var_col, i) for i in range(1,10)]
            for f in df_tsa[file_col]:
                print 'Extracting for array %s of approximately %s from:\n%s\n'\
                % (c, n_files, f)
                # If extract from the dataset depending on how year is stored
                dfs.append(extract_by_rowcol(df_xy, f, file_col, var_col, data_band,
                                             mosaic_tx, val_cols, data_type))
                c += 1
            # Comnbine all the pieces for this year
            df_yr = pd.concat(dfs)
            df_var[df_yr.columns] = df_yr
        
        # Write df_var with all years for this predictor to txt file
        this_bn = '%s_%s_kernelstats.txt' % (xy_txt_bn[:-4], var_name)
        this_txt = os.path.join(out_dir, this_bn)
        df_var.to_csv(this_txt, sep='\t')
        
    mosaic_ds = None
    
    # If out_text is specified, write the dataframe to a text file
    if out_dir:
        file_cols = [col for col in df_xy.columns if 'file' in col]
        out_bn = xy_txt_bn.replace('.txt', '_predictors.txt')
        out_txt = os.path.join(out_dir, out_bn)
        df_out = df_xy.drop(['tsa_id'] + file_cols, axis=1)# May want to keep row/col
        df_out.to_csv(out_txt, sep='\t')
        print 'Dataframe written to:\n', out_txt
    
    #return df_xy, df_tsa'''