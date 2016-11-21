# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 16:27:04 2016

@author: shooper
"""

import os
import sys
import ogr
import re
from glob import glob


def search_metadata(path, keyword, sep=':'):
    
    
    path = glob(path)[0]
    with open(path) as f:
        lines = [line.strip(' ').strip('\n') for line in f if keyword in line]
    
    val = lines[0].split(sep)[1]

    return val


def find_numbers(string):
    
    matches = re.findall('\d+\.*\d*', string)
    
    if len(matches) == 0:
        print 'No numbers found in %s. Returning 0...' % string
        return 0.0
    
    else:
        return [float(n) for n in matches]


def get_accuracy(string):
    
    numbers = find_numbers(string)
    print numbers
    
    if numbers == 0:
        return numbers
    
    elif len(numbers) == 1:            
        accuracy = numbers[0]
        if accuracy < 1:
            # It's a percent as decimal
            accuracy *= 100
        if 'error' in string.lower() and accuracy != 0:
            accuracy = 100 - accuracy
    
    # If more than 1 number was given take the mean
    else:
        accuracy = sum(numbers)/len(numbers)
    
    return accuracy


def main(search_str, shp, zone_field, keyword):
    
    ds = ogr.Open(shp, 1)
    lyr = ds.GetLayer()
    lyr_def = lyr.GetLayerDefn()
    
    features = [lyr.GetFeature(i) for i in range(lyr.GetFeatureCount())]
    zones = [(feat.GetFID(), feat.GetField(zone_field)) for feat in features]
    files = [(fid, search_str.format(z)) for fid, z in zones]
    acc_field = 'cnpy_acc'
    meta_field = 'cnpymeta'
    field_names = [lyr_def.GetFieldDefn(i).GetName() for i in range(lyr_def.GetFieldCount())]
    if not acc_field in field_names:
        acc_def = ogr.FieldDefn(acc_field, ogr.OFTReal)
        acc_def.SetPrecision(1)
        lyr.CreateField(acc_def)
    if not meta_field in field_names:
        meta_def = ogr.FieldDefn(meta_field, ogr.OFTString)
        meta_def.SetWidth(50)
        lyr.CreateField(meta_def)
    
    for fid, f in files:
        string = search_metadata(f, keyword)
        accuracy = get_accuracy(string)
        feat = lyr.GetFeature(fid)
        feat.SetField(acc_field, accuracy)
        feat.SetField(meta_field, string.strip('\r'))
        lyr.SetFeature(feat)
        feat.Destroy()
    
    ds.Destroy()
    

search_str = '/vol/v2/stem/canopy/truth_map/nlcd2001_canopy_mosaic_1-29-08/nlcd2001_canopy_metadata/metadata/zone{0}_canopy_meta_web_final/z*{0}*metadata_canopy.txt'
shp = '/vol/v2/stem/extent_shp/NLCD_modeling_regions.shp'
zone_field = 'Zone_ID'
keyword = 'Attribute_Accuracy_Value'

main(search_str, shp, zone_field, keyword)