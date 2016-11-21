# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 14:19:58 2015

@author: shooper
"""


def replace_val_with_array(tsa_ar, data_ar, tsa_id, offset):
    """
    Replace tsa_ar in place where tsa_id is replaced with data_ar where
    the two overlap according to offset
    """
    xoffset, yoffset = offset
    
    # If the upper left of the kernel is left or above the ul pixel, respectively
    if xoffset < 0 or yoffset < 0:
        data_cols = data_ar.shape[1]
        data_rows = data_ar.shape[0]    
        
        # Adjust the rows and columns to extract from the data only overlapping pixels 
        if xoffset < 0: 
            data_cols += xoffset # Subtract offset from last col to get from ds
            #xoffset = 0 # No offset because the array now starts at [0,0]
        if yoffset < 0: 
            data_rows += yoffset
            #yoffset = 0
            
        # Get just the portions of the data array that overlap
        data_ar = data_ar[abs(xoffset) :, abs(yoffset) :]
        
        # Replace tsa_id with the data array
        tsa_view = tsa_ar[: data_rows + 1, : data_cols + 1]
        tsa_mask = tsa_view == tsa_id
        np.copyto(tsa_view, data_ar, where=tsa_mask)
    
    # If the lower right of the kernel is to to the right or below the lr pixel
    elif xoffset + data_ar > xsize or yoffset + rows > ysize:
        data_cols = data_ar.shape[1]
        data_rows = data_ar.shape[0]    
        
        if xoffset + cols > xsize:
            data_cols -= xoffset 
        if yoffset + rows >  ysize:
            data_rows -= yoffset
     
        # Get just the portions of the data array that overlap
        data_ar = data_ar[: data_rows + 1, : data_cols + 1]
        
        # Replace tsa_id with the data array
        tsa_view = tsa_ar[data_rows :, data_cols :]
        tsa_mask = tsa_view == tsa_id
        np.copyto(tsa_view, data_ar, where=tsa_mask)

    # If they're miraculously the same size and in the same place
    else: np.copyto(tsa_ar, data_ar, where=(tsa_ar == tsa_id))
         
    return tsa_ar
    
    ''' Check if there's only 1 tsa_id and if so, get that dataset only
    and extract a kernel'''
    '''ul_x = x - (cols/2 * mosaic_tx[1])
    ul_y = y - (rows/2 * mosaic_tx[5])
    tx = ul_x, mosaic_tx[1], 0.0, ul_y, 0.0, mosaic_tx[5]
    print tx
    array_to_raster(tsa_ar, tx, prj, driver, out_path.replace('.bsq', '_tsa.bsq'), GDT_Int32)
    
    l = sorted(tsa_ids)
    tsa_ar[tsa_ar==l[0]] = 0
    tsa_ar[tsa_ar==l[1]] = 255/2
    tsa_ar[tsa_ar==l[2]] = 255
    
    cv2.imwrite(out_path.replace('.bsq', '.png'), tsa_ar)'''
    
def calc_offset(mosaic_tx, data_tx, tsa_ar_offset):
    '''
    Return the x and y offsets for a data and a TSA array, both referenced to mosaic_tx
    '''
    # Calc offset from global mosaic
    data_xoff = int((data_tx[0] - mosaic_tx[0])/mosaic_tx[1])
    data_yoff = int((data_tx[3] - mosaic_tx[3])/mosaic_tx[5])
    
    # Calculate the offset from the TSA array
    xoff = data_xoff - tsa_ar_offset[1]
    yoff = data_yoff - tsa_ar_offset[0]
    
    return pd.Series((xoff, yoff))

# ##From get_mosaic()
    prj = mosaic_ds.GetProjection()
    driver = mosaic_ds.GetDriver()

    '''df['dataset'] = df['file'].apply(lambda x: gdal.Open(x))
    df['data_tx'] = df['dataset'].apply(lambda x: x.GetGeoTransform())
    df[['xoff', 'yoff']] = df['data_tx'].apply(lambda x: calc_offset((ul_x, ul_y), x))
    
    # Get an array for each TSA
    df['data_array'] = df['dataset'].apply(
    lambda x: x.GetRasterBand(data_band).ReadAsArray())'''
    
    
def extract_kernel(ds, ar_coords, transform, xsize, ysize, nodata=0):
    # Modified original code from Zhiqiang Yang (read_spectral) at Oregon State University
    ''' 
    Return an array of pixel values from ds_ar centered around [x,y] of size (rows,cols)
    '''
    ul_x, ul_y, lr_x, lr_y = ar_coords
    #xoffset = int((x - transform[0])/transform[1]) - cols/2
    #yoffset = int((y - transform[3])/transform[5]) - rows/2 
    rows = int((lr_x - ul_x)/transform[1])
    cols = int((lr_y - ul_y)/transform[5])
    xoffset = int((ul_x - transform[0])/transform[1])
    yoffset = int((ul_y - transform[3])/transform[5])
    
    ds_cols = cols
    ds_rows = rows
    
    # Check if the array will be outside the bounds of the dataset
    x_dif = xoffset + cols - xsize
    y_dif = yoffset + rows - ysize
    if x_dif > 0: ds_cols -= x_dif
    if y_dif > 0: ds_rows -= y_dif

    # If the upper left of the kernel is left or above the ul pixel, respectively
    if xoffset < 0 or yoffset < 0:   
      
      # Adjust the rows and columns to extract only overlapping pixels from the data
        if xoffset < 0: 
            ds_cols = cols + xoffset # Subtract offset from last col to get from ds
            xoffset = 0 # No offset because the array now starts at [0,0]
        if yoffset < 0: 
            ds_rows = rows + yoffset
            yoffset = 0
        
        # Read the data as an array if necessary
        ds_ar = ds.ReadAsArray(xoffset, yoffset, ds_cols, ds_rows)
    
        # Fill an array of nodata with the extracted array to return an 
        #   array of the same size as the requested kernel
        ar = np.full((rows, cols), nodata, dtype=np.int)
        ar[rows - ds_rows:, cols - ds_cols:] = ds_ar
    
        return ar, (xoffset, yoffset)
    
    # If the lower right of the kernel is to to the right or below the lr pixel
    elif xoffset + cols > xsize or yoffset + rows > ysize:
    
        # Already adjusted the array size above so just extract the 
        #   array and do the same as above
        ds_ar = ds.ReadAsArray(xoffset, yoffset, ds_cols, ds_rows)
    
        ar = np.full((rows, cols), nodata, dtype=np.int)
        ar[: ds_rows, : ds_cols] = ds_ar
        
        return ar, (xoffset, yoffset)
    
    else: ar = ds.ReadAsArray(xoffset, yoffset, cols, rows) 
    
    return ar, (xoffset, yoffset)
    
    
def extract_kernel(array, ar_coords, transform, xsize, ysize, nodata=0):
    # Modified original code from Zhiqiang Yang (read_spectral) at Oregon State University
    ''' 
    Return an array of pixel values from ds_ar centered around [x,y] of size (rows,cols)
    '''
    ul_x, ul_y, lr_x, lr_y = ar_coords
    #xoffset = int((x - transform[0])/transform[1]) - cols/2
    #yoffset = int((y - transform[3])/transform[5]) - rows/2 
    cols = int((lr_x - ul_x)/transform[1])
    rows = int((lr_y - ul_y)/transform[5])
    xoffset = int((ul_x - transform[0])/transform[1])
    yoffset = int((ul_y - transform[3])/transform[5])
    
    ds_cols = cols
    ds_rows = rows
    
    # Check if the array will be partially outside the bounds of the dataset
    x_dif = xoffset + cols - xsize
    y_dif = yoffset + rows - ysize
    if x_dif > 0: ds_cols -= x_dif
    if y_dif > 0: ds_rows -= y_dif
        
    #Check if the array is entirely outside the bounds of the dataset
    #if xoffset + cols 
    

    # If the upper left of the kernel is left or above the ul pixel, respectively
    if xoffset < 0 or yoffset < 0:   
      
      # Adjust the rows and columns to extract only overlapping pixels from the data
        if xoffset < 0: 
            ds_cols = cols + xoffset # Subtract offset from last col to get from ds
            xoffset = 0 # No offset because the array now starts at [0,0]
        if yoffset < 0: 
            ds_rows = rows + yoffset
            yoffset = 0
        
        # Get the portion of the array that the kernel covers
        ds_ar = array[yoffset : yoffset + ds_rows,  xoffset : xoffset + ds_cols]
    
        # Fill an array of nodata with the extracted array to return an 
        #   array of the same size as the requested kernel
        ar = np.full((rows, cols), nodata, dtype=np.int)
        ar[rows - ds_rows:, cols - ds_cols:] = ds_ar
    
        return ar, (xoffset, yoffset)
    
    # If the lower right of the kernel is to to the right or below the lr pixel
    elif xoffset + cols > xsize or yoffset + rows > ysize:
    
        # Already adjusted the array size above so just extract the 
        #   array and do the same as above
        ds_ar = array[yoffset : yoffset + ds_rows,  xoffset : xoffset + ds_cols]
        #ds_ar = ds.ReadAsArray(xoffset, yoffset, ds_cols, ds_rows)
    
        ar = np.full((rows, cols), nodata, dtype=np.int)
        ar[: ds_rows, : ds_cols] = ds_ar
        
        return ar, (xoffset, yoffset)
    
    else: ar = array[yoffset : yoffset + rows,  xoffset : xoffset + cols] 
    
    return ar, (xoffset, yoffset)

''' This function doesn't work if the the ul and lr are beyond the boundaries of the dataset'''
def extract_kernel(ds, data_band, ar_coords, transform, xsize, ysize, nodata=0):
    # Modified original code from Zhiqiang Yang (read_spectral) at Oregon State University
    ''' 
    Return an array of pixel values from ds_ar centered around [x,y] of size (rows,cols)
    '''
    ul_x, ul_y, lr_x, lr_y = ar_coords
    xsize -= 1 #xsize and ysize are number of row/cols, not index of last row/col
    ysize -= 1
    #xoffset = int((x - transform[0])/transform[1]) - cols/2
    #yoffset = int((y - transform[3])/transform[5]) - rows/2 
    cols = int((lr_x - ul_x)/transform[1])
    rows = int((lr_y - ul_y)/transform[5])
    xoffset = int((ul_x - transform[0])/transform[1])
    yoffset = int((ul_y - transform[3])/transform[5])
    
    ds_cols = cols
    ds_rows = rows
    print xoffset
    print yoffset
    # Check if the array will be outside the bounds of the dataset
    x_dif = xoffset + cols - xsize 
    y_dif = yoffset + rows - ysize
    if x_dif > 0: ds_cols -= x_dif
    if y_dif > 0: ds_rows -= y_dif

    # If the upper left of the kernel is left or above the ul pixel, respectively
    if xoffset < 0 or yoffset < 0:   
      
      # Adjust the rows and columns to extract only overlapping pixels from the data
        if xoffset < 0: 
            ds_cols = cols + xoffset # Subtract offset from last col to get from ds
            xoffset = 0 # No offset because the array now starts at [0,0]
        if yoffset < 0: 
            ds_rows = rows + yoffset
            yoffset = 0
        
        # Read the data as an array if necessary
        ds_ar = ds.GetRasterBand(data_band).ReadAsArray(xoffset, yoffset, ds_cols, ds_rows)
    
        # Fill an array of nodata with the extracted array to return an 
        #   array of the same size as the requested kernel
        ar = np.full((rows, cols), nodata, dtype=np.int32)
        ar[rows - ds_rows:, cols - ds_cols:] = ds_ar
    
        return ar, (xoffset, yoffset)
    
    # If the lower right of the kernel is to to the right or below the lr pixel
    elif xoffset + cols > xsize or yoffset + rows > ysize:
    
        # Already adjusted the array size above so just extract the 
        #   array and do the same as above
        ds_ar = ds.GetRasterBand(data_band).ReadAsArray(xoffset, yoffset, ds_cols, ds_rows)
    
        ar = np.full((rows, cols), nodata, dtype=np.int32)
        ar[: ds_rows, : ds_cols] = ds_ar
        
        return ar, (xoffset, yoffset)
    
    else: ar = ds.GetRasterBand(data_band).ReadAsArray(xoffset, yoffset, cols, rows) 
    
    return ar, (xoffset, yoffset)