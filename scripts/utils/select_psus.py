# -*- coding: utf-8 -*-
"""
Select primary sampling units (PSUs) from CONUS tiles for evaluating CONUS maps

example command:
python select_psus.py /vol/v1/general_files/datasets/spatial_data/conus_tile_system/conus_tile_system_15_sub_epsg5070.shp /vol/v1/proj/stem_improv_paper/vector/NA_CEC_Eco_Level2_clipped_to_conus_epsg5070.shp 100 /vol/v2/stem/conus/vector/sampling_psu.shp

"""

import os, sys, random
import pandas as pd
import numpy as np
from osgeo import ogr

package_dir = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(package_dir)
from stem import coords_to_shp
from lthacks import attributes_to_df, createMetadata


def main(tile_shp, strata_shp, n_psu, out_shp, strata_id_field='NA_L2NAME', min_dist=173779, split_tags=['2001','2011']):
    
    n_psu = int(n_psu)
    
    tiles = attributes_to_df(tile_shp)
    tile_ds = ogr.Open(tile_shp)
    tile_lyr = tile_ds.GetLayer()
    tiles['xctr'] = (tiles.xmax -tiles.xmin)/2 + tiles.xmin
    tiles['yctr'] = (tiles.ymax -tiles.ymin)/2 + tiles.ymin
    tiles['ul_x'] = tiles.xmin
    tiles['lr_x'] = tiles.xmax
    tiles['ul_y'] = tiles.ymax
    tiles['lr_y'] = tiles.ymin
    
    strata = attributes_to_df(strata_shp)
    #
    strata_ds = ogr.Open(strata_shp)
    strata_lyr = strata_ds.GetLayer()
    
    # Get areas and calculate proportions of total
    for feat in strata_lyr:
        fid = feat.GetFID()
        geom = feat.GetGeometryRef()
        area = geom.GetArea()
        strata.loc[fid, 'area'] = area
    
    # Features could be multipart, so calculate sums for all parts of same stratum
    unique_names = strata[strata_id_field].unique()
    summed_areas = pd.Series({name: strata.loc[strata[strata_id_field] == name, 'area'].sum() for name in unique_names if name != 'WATER'})
    strata.drop_duplicates(strata_id_field, inplace=True)
    strata.set_index(strata_id_field, inplace=True)
    strata.drop('WATER', inplace=True)
    
    strata['area'] = summed_areas/summed_areas.sum()
    strata['n_psu'] = (strata.area * n_psu).round().astype(int)
    strata.loc[strata.n_psu == 0, 'n_psu'] = 1
    
    # Randomly shuffle strata so the same strata don't always influence availble
    #   psus
    strata = strata.sample(frac=1) 
    candidates = tiles.copy()
    fids = []
    strata_names = {}
    for i, (stratum_name, stratum) in enumerate(strata.iterrows()):
        print i, stratum_name, ':',
        strata_lyr.SetAttributeFilter("%s = '%s'" % (strata_id_field, stratum_name))
        strata_feat = strata_lyr.GetNextFeature()
        strata_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        while strata_feat:
            g = strata_feat.GetGeometryRef()
            strata_geom = strata_geom.Union(g)
            strata_feat = strata_lyr.GetNextFeature()
          
        # find all tile features that intersect this stratum
        overlapping = []
        print 'getting overlapping...',
        for t_fid in candidates.index:
            tile_feature = tile_lyr.GetFeature(t_fid)
            tile_geom = tile_feature.GetGeometryRef()
            if strata_geom.Intersects(tile_geom):
                overlapping.append(t_fid)#'''
        if len(overlapping) == 0:
            continue
        
        print 'selecting...\n'
        for j in range(stratum.n_psu):
            this_fid = random.sample(overlapping, 1)
            fids.extend(this_fid)
            selected = tiles.loc[fids]
            strata_names[this_fid[0]] = stratum_name
            for ti, c_tile in candidates.iterrows():   
                if np.any(np.sqrt((selected.xctr - c_tile.xctr)**2 + (selected.yctr - c_tile.yctr)**2) <= min_dist):
                    candidates.drop(ti, inplace=True)
                    # Additionally remove tiles from overlapping list so they're not selected
                    if ti in overlapping:
                        # Might not be depending on search distance
                        overlapping.remove(ti)
    
    selected[strata_id_field] = pd.Series(strata_names)
    
    if split_tags:
        #random_ids = random.sample(selected.index, strata.n_psu.sum()/2)
        selected1 = selected.sample(frac=.5)
        selected2 = selected.loc[~selected.index.isin(selected1.index)]
        coords_to_shp(selected1, tile_shp, out_shp.replace('.shp', '_%s.shp') % split_tags[0])
        coords_to_shp(selected2, tile_shp, out_shp.replace('.shp', '_%s.shp') % split_tags[1])
    else:
        #selected.to_csv(out_shp.replace('.shp', '.txt'))
        coords_to_shp(selected, tile_shp, out_shp)
    
    strata_ds, strata_lyr, strata_feat = None, None, None
    tile_ds, tile_lyr = None, None

if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
            
                
                                