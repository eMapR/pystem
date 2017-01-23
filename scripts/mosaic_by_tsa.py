
import os
import sys
import fnmatch
import pandas as pd
import numpy as np
from osgeo import gdal
from gdalconst import *


def read_params(params):
    '''
    Read parameter file and parse into dictionary
    '''
    if not os.path.exists(params):
        print 'Parameter file given does not exist:\n%s' % params
        return None
    
    d = {}
    try:
        with open(params) as f:
            input_vars = [line.split(";") for line in f]
    except: 
        print 'Problem reading parameter file:\n%s' % params
        return None
    
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
            '"%s"' % var[1].strip(" ").replace("\n", "")

    print 'Parameters read from:\n%s\n' % params
    
    return d
    
    
def find_file(basepath, tsa_str, search_str, path_filter=None):
    '''
    Return the full path within the directory tree /baspath/tsa_str if search_str
    is in the filename. Optionally, if path_filter is specified, only a path that
    contains path_filter will be returned.
    '''
    """if not os.path.exists(basepath):
        print 'basepath does not exist: \n%s' % basepath
        return None
     
    bp = os.path.join(basepath, tsa_str)
    
    # Search the bp directory tree. If search_str is in a file, get the full path.
    paths = []
    for root, dirs, files in os.walk(bp, followlinks=True):
        these_files = [f for f in files if fnmatch.fnmatch(f, search_str)]
        if len(these_files) > 0:
            these_paths = [os.path.join(root, f) for f in these_files]
            paths.extend(these_paths)
    
    if path_filter:
        paths = [p for p in paths if path_filter in p]
    
    '''if len(paths) > 1:
        print 'Multiple files found for tsa: %s' % tsa_str
        for p in paths:
            print p
        return None'''

    return paths[0]"""

    if not os.path.exists(basepath):
        print 'basepath does not exist: \n%s' % basepath
        return None
     
    bp = os.path.join(basepath, tsa_str)

    # Search the bp directory tree. If search_str is in a file, get the full path.
    paths = []
    for root, dirs, files in os.walk(bp, followlinks=True):
        these_paths = [os.path.join(root, f) for f in files]
        these_paths = fnmatch.filter(these_paths, search_str)
        if len(these_paths) > 0:
            paths.extend(these_paths)
    
    # If path filter was specified, remove any paths that contain it
    '''if not path_filter == '':
        [paths.remove(p) for p in paths if fnmatch.fnmatch(p, path_filter)]
        #paths = [p for p in paths if fnmatch.fnmatch(p, path_filter)]'''
    
    if len(paths) > 1:
        '''print 'Multiple files found for tsa: ' + tsa_str
        for p in paths:
            print p
        print 'Selecting the first one found...\n' '''
        #return None
        
    if len(paths) < 1:
        sys.exit(('No files found for tsa {0} with basepath {1} and ' +\
        'search_str {2}\n').format(tsa_str, basepath, search_str))
        #return None
        
    return paths[0]


def calc_offset(tsa_ul, data_tx):
    '''
    Return the row and col offset of a data array from a tsa_array
    '''
    tsa_x, tsa_y = tsa_ul
    data_x, data_y = data_tx[0], data_tx[3]
    
    row_off = int((data_y - tsa_y)/data_tx[5])
    col_off = int((data_x - tsa_x)/data_tx[1])
    
    #return pd.Series((row_off, col_off))
    return row_off, col_off


def get_array(filepath, data_band, coords):
    '''
    Return an array, geo transform, and offset from ul_off from the 
    raster dataset with the path, filepath
    '''

    ds = gdal.Open(filepath)
    tx = ds.GetGeoTransform()
    
    if len(coords) == 2:
        ar = ds.GetRasterBand(data_band).ReadAsArray()
        row_off, col_off = calc_offset(coords, tx)
    
    else:
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        ar, offsets = extract_kernel(ds, data_band, coords, tx, xsize, ysize)
        row_off, col_off = offsets
        
    return tx, ar, row_off, col_off


def get_offset_array_indices(ar1_shape, ar2_shape, ar2_offset):
    '''
    Return the indices to broadcast ar2 into ar1
    '''
    # Initialize stuff
    yoffset, xoffset = ar2_offset
    ar2_rows, ar2_cols = ar2_shape
    ar1_rows, ar1_cols = ar1_shape
    # Initialize indices for arrays in case they don't need to be set by
    #  the offset
    ar2_col_l, ar2_row_u, ar1_col_l, ar1_row_u = 0, 0, 0, 0
    ar2_col_r, ar2_row_d = ar2_cols, ar2_rows
    ar1_col_r, ar1_row_d = ar1_cols, ar1_rows
    
    # If the upper left pixel of the ar2 is left of the ul pixel 
    #   of ar1, reset the index of the ar1 to get only cols to the right
    #   of the offset
    if xoffset < 0:
        ar2_col_l = abs(xoffset)
    # Otherwise, do the same for ar1
    else:
        ar1_col_l = xoffset
        
    # Do the same checks and adjustments if the ar2 ul is above the ar1 ul
    if yoffset < 0:
        ar2_row_u = abs(yoffset)
    else:
        ar1_row_u = yoffset
    
    # If the lower right pixel of ar2 is to the right of the lr pixel of
    #   ar1, truncate ar2 by the difference
    if xoffset + ar2_cols > ar1_cols:
        ar2_col_r = ar1_cols - xoffset
        #ar2_col_r = -(xoffset + ar2_cols - ar1_cols)
    # Otherwise, the ar1 array needs to be truncated
    else:
        ar1_col_r = ar2_cols + xoffset

    # Do the same if either lr is below the other
    if yoffset + ar2_rows > ar1_rows:
        ar2_row_d = ar1_rows - yoffset
        #ar2_row_d = -(yoffset + ar2_rows - ar1_rows)
    # Otherwise, the ar1 array needs to be truncated
    else:
        ar1_row_d = ar2_rows + yoffset

    ar1_bounds = ar1_row_u, ar1_row_d, ar1_col_l, ar1_col_r
    ar2_bounds = ar2_row_u, ar2_row_d, ar2_col_l, ar2_col_r
    
    return ar1_bounds, ar2_bounds


def extract_kernel(ds, data_band, ar_coords, transform, xsize, ysize, nodata=0):
    
    ul_x, ul_y, lr_x, lr_y = ar_coords #Projected coords
    #xsize -= 1 #xsize and ysize are number of row/cols, not index of last row/col
    #ysize -= 1

    cols = int((lr_x - ul_x)/transform[1])
    rows = int((lr_y - ul_y)/transform[5])
    
    xoffset = int((ul_x - transform[0])/transform[1])
    yoffset = int((ul_y - transform[3])/transform[5])
    
    # Get offsets and number of rows/cols for the data array
    ds_inds, array_inds = get_offset_array_indices((ysize, xsize), (rows, cols), (yoffset, xoffset))
    ds_row_ul, ds_row_lr, ds_col_ul, ds_col_lr = ds_inds
    ds_cols = ds_col_lr - ds_col_ul
    ds_rows = ds_row_lr - ds_row_ul
    band = ds.GetRasterBand(data_band)
    ar_ds = band.ReadAsArray(ds_col_ul, ds_row_ul, ds_cols, ds_rows)
    
    # Broadcast ds_ar into an array of nodata values
    ar = np.full((rows, cols), nodata, dtype=np.int32)
    ar_row_ul, ar_row_lr, ar_col_ul, ar_col_lr = array_inds
    ar[ar_row_ul : ar_row_lr, ar_col_ul : ar_col_lr] = ar_ds
    
    return ar, (ar_row_ul, ar_col_ul)
    
    
def replace_val_with_array(tsa_ar, data_ar, tsa_id, offset):
    '''
    Replace tsa_ar in place where tsa_ar == tsa_id and where the two 
    overlap according to offset
    '''
    tsa_shape  = tsa_ar.shape
    data_shape = data_ar.shape
    tsa_bounds, data_bounds = get_offset_array_indices(tsa_shape, data_shape, offset)
    tsa_row_u,  tsa_row_d,  tsa_col_l,  tsa_col_r = tsa_bounds
    data_row_u, data_row_d, data_col_l, data_col_r = data_bounds
    
    # Replace tsa_id with the data array
    data_view = data_ar[data_row_u : data_row_d, data_col_l : data_col_r]
    tsa_view = tsa_ar[tsa_row_u : tsa_row_d, tsa_col_l : tsa_col_r]
    tsa_mask = tsa_view == tsa_id
    #print tsa_id
    #import pdb; pdb.set_trace()    
    np.copyto(tsa_view, data_view, where=tsa_mask)
         

def array_to_raster(array, tx, prj, driver, out_path, dtype, nodata=None):
    # From intersectMask.py
    '''
    Save a numpy array as a new raster
    '''
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
    

#def get_mosaic(mosaic_tx, tsa_ar, tsa_off, ar_coords, tsa_txt,  data_band, files=None, basepath=None, search_str=None, path_filter=None, out_path=None, prj=None, driver=None):
def get_mosaic(mosaic_tx, tsa_ids, tsa_ar, ar_coords,  data_band, files=None, **kwargs):

    ul_x, ul_y, lr_x, lr_y = ar_coords

    tsas = zip(tsa_ids, [int(tsa) for tsa in tsa_ids])
    df = pd.DataFrame(tsas, columns=['tsa_str','tsa_id'])
    #import pdb; pdb.set_trace()
    
    # If files weren't given, store the filepath for each TSA in the df
    if not files:
        df['file'] = df['tsa_str'].apply(
        lambda x: find_file(basepath, x, search_str, path_filter))
        if df['file'].isnull().any():
            print 'Could not find unique file for TSA. Try a different ' +\
            'search string or specify a path filter'
            return None
        these_files = df.file
    else:
        these_files = [f for f in files for tsa in df['tsa_str'] if tsa in f]
    
 
    # Store an array for each TSA and calculate offsets from tsa_ar    
    data_ars = [get_array(f, data_band, (ul_x, ul_y)) for f in these_files]
    df_ars = pd.DataFrame(data_ars, index=df.index, 
                          columns=['data_tx', 'data_array', 'yoff', 'xoff'])
    df[['data_tx', 'data_array', 'yoff', 'xoff']] = df_ars
    
    # Replace each tsa_id in tsa_ar with data values. But first, make 
    #   tsa_ar of the same dtype as the data arrays.
    #tsa_ar = tsa_ar.astype(df.ix[df.index.min(), 'data_array'].dtype)
    
    [replace_val_with_array(tsa_ar, row['data_array'], row['tsa_id'],
                            (row['yoff'],row['xoff'])) 
                            for ind, row in df.iterrows()]
    
    # Save the new mosaic as a raster if out_path is specified
    if 'out_path' in locals():
        tx = ul_x, mosaic_tx[1], 0.0, ul_y, 0.0, mosaic_tx[5]
        array_to_raster(tsa_ar, tx, prj, driver, out_path, GDT_Int32) 
        ''' fix dtype so that it can be a float if the array is'''
    del df_ars, df
    
    return tsa_ar.astype(np.int32)


def main(mosaic_path, basepath, search_str, ar_coords, data_band):
    
    # Read in the lookup table of TSAs
    #df = pd.read_csv(tsa_txt, sep='\t', dtype={'tsa_str':object})
    #df = df[df['tsa_id'].isin(tsa_ids)]

    ds_tsa = gdal.Open(mosaic_path)
    tx_tsa = ds_tsa.GetGeoTransform()
    xsize = ds_tsa.RasterXSize - 1
    ysize = ds_tsa.RasterYSize - 1
    
    ar_tsa, tsa_off = extract_kernel(ds_tsa, 1, ar_coords, tx_tsa, xsize, ysize)
    ds_tsa = None
    tsa_ids = np.unique(ar_tsa)
    tsa_strs = ['0' + str(tsa) for tsa in tsa_ids if tsa!=0]

    ar = get_mosaic(tx_tsa, tsa_strs, ar_tsa, ar_coords, data_band)
    
    return ar

    
''' ###### Testing ###### '''


##x, y = -1900000, 2770000
#x, y = -2040000, 2750000
#x2, y2 = -1830000, 2675000
#cols = 4000
#rows = 4000
#tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'


"""mosaic_path = '/vol/v1/general_files/datasets/spatial_data/CAORWA_TSA_lt_only.bsq'
'''ul_x = -1628643
ul_y = 3203960
lr_x = -1328643
lr_y = 2803960 '''                                

basepath = '/vol/v1/scenes'
#basepath = '/vol/v3/scenes'
search_str = '*_mse_split.bsq'
#search_str = '*2010_*ledaps.bsq'
path_filter = None
out_dir = '/home/server/student/homes/shooper'
data_band = 18

ds = gdal.Open(mosaic_path)
ar_tsa = ds.ReadAsArray()
tx = ds.GetGeoTransform()
xsize = ds.RasterXSize
ysize = ds.RasterYSize
ul_x, x_res, _, ul_y, _, y_res = tx
ar_coords = [ul_x,
             ul_y,
             ul_x + (xsize * x_res), 
ul_y + (ysize * y_res)]
#ar, off = extract_kernel(ds, 1, ar_coords, tx, xsize, ysize)

prj = ds.GetProjection()
driver = ds.GetDriver()
m_ulx, x_res, x_rot, m_uly, y_rot, y_res = tx
tx1 = ul_x, x_res, x_rot, ul_y, y_rot, y_res
ds = None
#out_dir = '/vol/v2/stem/scripts/testing'
out_path = os.path.join(out_dir, 'mse_CAORWA.bsq')

ar = main(mosaic_path, basepath, search_str, ar_coords, data_band)
nodata_mask = ar_tsa == 0
ar[nodata_mask] == -9999
array_to_raster(ar, tx, prj, driver, out_path, GDT_Int16, nodata=-9999)#"""

    
    
    
    

    
    
    
