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

def kernel_from_shp(mosaic_lyr, coords, mosaic_tx, nodata, val_field='id'):
    
    # Get rid of the annoying warning about spatial refs when rasterizing
    #gdal.PushErrorHandler('CPLQuietErrorHandler')
    
    ul_x, x_res, _, ul_y, _, y_res = mosaic_tx

    #set_wkt = 'POLYGON (({0} {1}, {2} {3}, {4} {5}, {6} {7}))'.format(coords.ul_x,coords.ul_y, coords.lr_x, coords.ul_y, coords.lr_x, coords.lr_y, coords.ul_x, coords.lr_y)
    support_set = Polygon([[coords.ul_x, coords.ul_y], 
                                           [coords.lr_x, coords.ul_y],
                                           [coords.lr_x, coords.lr_y],
                                           [coords.ul_x, coords.lr_y]
                                           ])
    #support_set = wkt.loads(set_wkt)
    #support_set = ogr.CreateGeometryFromWkt(set_wkt)
    #support_set.CloseRings()
    
    ''' figure out a better way to compute intersection than looping though'''
    #intersection = ogr.Geometry(ogr.wkbMultiPolygon)
    
    # Create the raster datasource and the intersection layer in memory
    srs = mosaic_lyr.GetSpatialRef()
    mem_driver = ogr.GetDriverByName('Memory')
    shp_driver = ogr.GetDriverByName('Esri Shapefile')
    ds_set = mem_driver.CreateDataSource('mem')
    ds_inter = mem_driver.CreateDataSource('inter')
    intersection_lyr = ds_inter.CreateLayer('intersection', None, ogr.wkbPolygon)
    '''set_lyr = ds_set.CreateLayer('set', srs, ogr.wkbPolygon)
    #tile_lyr = ds_mem.CreateLayer('tile', srs, ogr.wkbPolygon)
    lyr_def = mosaic_lyr.GetLayerDefn()
    set_feature = ogr.Feature(lyr_def)
    set_feature.SetGeometry(support_set)
    set_lyr.CreateFeature(set_feature.Clone())
    mosaic_lyr.Intersection(set_lyr, intersection_lyr)'''
    
    #out_shp = shp_driver.CreateDataSource('/vol/v2/stem/conus_testing/models/landcover_20170419_1443/intersection_delete.shp')
    #out_lyr = out_shp.CreateLayer('intersection_delete', mosaic_lyr.GetSpatialRef(), ogr.wkbMultiPolygon)
    
    lyr_def = mosaic_lyr.GetLayerDefn()
    for i in range(lyr_def.GetFieldCount()):
        field_def = lyr_def.GetFieldDefn(i)
        intersection_lyr.CreateField(field_def)
        #out_lyr.CreateField(field_def)
    #for i in range(mosaic_lyr.GetFeatureCount()):
    for feature in mosaic_lyr:
        #feature = mosaic_lyr.GetFeature(i)
        #geometry = feature.GetGeometryRef()
        geometry = WKT.loads(feature.GetGeometryRef().ExportToWkt())
        print geometry.is_valid
        #geometry.CloseRings()
        '''if geometry.Intersects(support_set):
            this_intersection = ogr.Geometry(ogr.wkbMultiPolygon)
            for i in range(geometry.GetGeometryCount()):
                g = geometry.GetGeometryRef(i)
                if g.Intersects(support_set):
                    intersected_geom = g.Intersection(support_set)
                    #import pdb; pdb.set_trace()
                    try:
                        this_intersection.AddGeometry(intersected_geom)
                    except:
                        import pdb; pdb.set_trace()#'''
        if geometry.intersects(support_set):
            #intersected_geom = geometry.Intersection(geometry)
            intersected_wkt = geometry.intersection(support_set).wkt
            intersected_geom = ogr.CreateGeometryFromWkt(intersected_wkt)
            intersected_feature = ogr.Feature(lyr_def)
            intersected_feature.SetGeometry(intersected_geom)
            #import pdb; pdb.set_trace()
            for i in range(lyr_def.GetFieldCount()):
                name = lyr_def.GetFieldDefn(i).GetName()
                val = feature.GetField(i)
                intersected_feature.SetField(name, val)
            #out_lyr.CreateFeature(intersected_feature.Clone())
            intersection_lyr.CreateFeature(intersected_feature.Clone())
            feature.Destroy()
            intersected_feature.Destroy()
    #out_shp.Destroy()#"""
    
    
    #Get number of rows and cols of the output array and the input layer
    x1, x2, y1, y2 = support_set.GetEnvelope()
    cols = abs(int((x2 - x1)/x_res))
    rows = abs(int((y2 - y1)/y_res))
    
    m_xmin, m_xmax, m_ymin, m_ymax = mosaic_lyr.GetExtent()
    xsize = abs(int((m_xmax - m_xmin)/x_res))
    ysize = abs(int((m_ymax - m_ymin)/y_res))
    
    tx = x1, x_res, 0, y2, 0, y_res #Assumes x > to the east, y > north
    offset = calc_offset((m_xmin, m_ymax), tx)
    
    ''' Figuer out how to burn nodata values directly with Rasterize() rather than with mask'''
    ras_driver = gdal.GetDriverByName('MEM')
    ds_ras = ras_driver.Create('', cols, rows, 1, gdal.GDT_Int16)
    ds_mask = ras_driver.Create('', cols, rows, 1, gdal.GDT_Byte)
    ds_ras.SetGeoTransform(tx)
    ds_mask.SetGeoTransform(tx)
    gdal.RasterizeLayer(ds_ras, [1], intersection_lyr, options=["ATTRIBUTE=%s" % val_field])
    gdal.RasterizeLayer(ds_mask, [1], intersection_lyr, burn_values=[1])
    ar = ds_ras.ReadAsArray()
    mask = ds_mask.ReadAsArray()
    import pdb; pdb.set_trace()
    mask = mask.astype(bool)
    ar[~mask] = nodata
    ds_ras = None
    ds_mask = None
    ds_mem.Destroy()

    _, array_inds = get_offset_array_indices((ysize, xsize), (rows, cols), offset)
    ul_row, ar_row_lr, ar_col_ul, ar_col_lr = array_inds
    offset = array_inds[2], array_inds[0]

    return ar, offset