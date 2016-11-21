# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 20:12:05 2016

@author: shooper
"""

def get_overlapping_polys(src_shp, ovr_shp, out_shp):
    '''
    Return a shapefile of all features in ds_src that touch ds_ovr
    '''
    ds_src = ogr.Open(src_shp)
    if ds_src == None:
        print 'Shapefile does not exist or is not valid:\n%s' % src_shp
        return None
    lyr_src = ds_src.GetLayer()
    srs_src = lyr_src.GetSpatialRef()
    
    ds_ovr = ogr.Open(ovr_shp)
    if ds_ovr == None:
        print 'Shapefile does not exist or is not valid:\n%s' % ovr_shp
        return None
    lyr_ovr = ds_ovr.GetLayer()
    
    # Create the output dataset
    driver = ds_src.GetDriver()
    if os.path.exists(out_shp):
        os.remove(out_shp)   
    ds_out = driver.CreateDataSource(out_shp)
    lyr_out = ds_out.CreateLayer(os.path.basename(out_shp)[:-4], srs_src, geom_type=ogr.wkbMultiPolygon) 
    lyr_out_def = lyr_out.GetLayerDefn()
    
    #Get field definitions
    lyr_src_def = lyr_src.GetLayerDefn()
    for i in range(lyr_src_def.GetFieldCount()):
        field_def = lyr_src_def.GetFieldDefn(i)
        lyr_out.CreateField(field_def) 
    
    # Get all of the overlap geometries up front
    ovr_geoms = []
    for i in xrange(lyr_ovr.GetFeatureCount()):
        feat_ovr = lyr_ovr.GetFeature(i)
        ovr_geoms.append(feat_ovr.GetGeometryRef())
        feat_ovr.Destroy()
    
    t0 = time.time()
    # Loop through each feautre and check for overlap
    feat_src = lyr_src.GetNextFeature()
    while feat_src:
        geom_src = feat_src.GetGeometryRef()
        
        for geom_ovr in ovr_geoms:
            #feat_ovr = lyr_ovr.GetFeature(j)
            #geom_ovr = feat_ovr.GetGeometryRef()
            
            # If there's any overlap, add the feature to the lyr_out
            if geom_ovr.Intersect(geom_src):
                feat_out = ogr.Feature(lyr_out_def)
                feat_out.SetGeometry(geom_src)
                # Get the fields from the source
                for j in range(lyr_out_def.GetFieldCount()):
                    feat_out.SetField(lyr_out_def.GetFieldDefn(i).GetNameRef(), feat_src.GetField(i))
                lyr_out.CreateFeature(feat_out)
                feat_out.Destroy()
                break
        
        feat_src.Destroy()
        feat_src = lyr_src.GetNextFeature()
    
    print 'Total time: %.1f\n' % (time.time() - t0)
    
    ds_src.Destroy()
    ds_ovr.Destroy()
    lyr_src = None
    lyr_ovr = None
    
    print 'Shapefile written to: ', out_shp  