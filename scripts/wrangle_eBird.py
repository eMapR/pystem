import os
import time
import subprocess
from osgeo import gdal, ogr, osr
import pandas as pd
import numpy as np

gdal.UseExceptions()


def parse_erd(txt, states, state_field, species, columns, chunk_size, nan_cols, dup_col='PRIMARY_CHECKLIST_FLAG',
              sep=',', na_vals=['?', 'NA'], no_count='X'):
    print txt
    if not os.path.exists(txt):
        print 'Input txt file does not exist:\n%s' % txt
        return False

    columns[species] = object  # Column is mostly ints but has some strings

    # Read in chunks and only grab the records from states in the states list
    reader = pd.read_csv(txt, usecols=columns.keys(), sep=sep, chunksize=chunk_size, na_values=na_vals, dtype=columns)
    df = pd.concat([chunk[chunk[state_field].isin(states)] for chunk in reader], ignore_index=True)
    df.ix[df[species] == no_count, species] = 1  # If no count was given, replace with a 1
    df[species] = df[species].astype(int)

    # Turn counts into 0-1 presence-absence
    genus, species_name = species.split('_')[0].lower()
    target_col = '_'.join(genus[0], species_name) # Response variable name
    df.ix[df[species] > 0, target_col] = 1
    df.ix[df[species] == 0, target_col] = 0
    #df.drop([species], axis=1, inplace=True)

    # Get rid of any records with nan values in necessary columns and duplicates
    #df = df[df[nan_cols].notnull().all(axis=1) & df[dup_col] == 1]
    return df


def reproject_coords(df, lat_field, lon_field, src_shp, dst_shp, out_dir):
    # Get spatial ref from source and destination shps
    src_ds = ogr.Open(src_shp)
    src_lyr = src_ds.GetLayer()
    src_srs = src_lyr.GetSpatialRef()
    src_ds, src_lyr = None, None
    src_wkt = src_srs.ExportToWkt()

    dst_ds = ogr.Open(dst_shp)
    dst_lyr = dst_ds.GetLayer()
    dst_srs = dst_lyr.GetSpatialRef()
    dst_ds, dst_lyr = None, None
    dst_wkt = dst_srs.ExportToWkt()

    # Write the input coords to a text files
    in_temp_txt = os.path.join(out_dir, 'latlon.txt')
    df[[lon_field, lat_field]].to_csv(in_temp_txt, sep=' ', index=False, header=False)

    # Reproject the input coords and write the result to a text file
    out_temp_txt = os.path.join(out_dir, 'xy_prj.txt')
    cmd = 'gdaltransform -s_srs {0} -t_srs {1} -output_xy < {2} > {3}'.format(src_wkt, dst_wkt, in_temp_txt,
                                                                              out_temp_txt)
    subprocess.call(cmd, universal_newlines=True, shell=True)
    os.remove(in_temp_txt)

    # Read the projected coords back in as a text file and merge with df
    df_xy = pd.read_csv(out_temp_txt, sep=' ', header=None, names=['x', 'y'])
    df[['x', 'y']] = df_xy
    df[['x', 'y']] = df[['x', 'y']].astype(int)  # Convert to int

    return df, out_temp_txt

def main(years, txt_template, states, species, columns, state_field, lat_field, lon_field, src_shp, dst_shp,
         out_dir, xy_txt=None, nan_cols=None, dup_col=None, drop_cols=None):
    
    t0 = time.time()
    if not os.path.exists(out_dir):
        print '\nout_dir does not currently exist:\n%s\nCreating output directory...' % out_dir
        os.mkdir(out_dir)

    t1 = time.time()
    print '\nConcatenating records from multiple years...'
    df = pd.concat([parse_erd(txt_template % yr, states, state_field, species, columns, 5000, nan_cols) for yr in years])
    print '%.1f seconds\n' % ((time.time() - t1)/60)

    if not xy_txt:
        t1 = time.time()
        print 'No projected coords provided. Calculating projected coords from lat lon...'
        df_xy, xy_txt = reproject_coords(df, lat_field, lon_field, src_shp, dst_shp, out_dir)
        print '%.1f seconds\n' % (time.time() - t1)
    else:
        try:
            df_xy = pd.read_csv(xy_txt)
            df[['x', 'y']] = df_xy
        except Exception as e:
            print e
            import pdb; pdb.set_trace()

    out_txt = os.path.join(out_dir, 'erd_v5_%s.txt' % species)
    #df.drop([lat_field, lon_field], axis=1, inplace=True)
    if nan_cols and dup_col:
        print 'Dropping duplicates and records with null values in \n%s\n' % '\n'.join(nan_cols)
        df = df[df[nan_cols].notnull().all(axis=1) & df[dup_col] == 1]
    columns = df.columns
    if drop_cols:
        drop_cols = [col for col in drop_cols if col in columns]
        import pdb; pdb.set_trace()
        df.drop(drop_cols, axis=1, inplace=True)
            
    df['obs_id'] = range(len(df))
    df.set_index('obs_id', inplace=True)
    df.to_csv(out_txt, sep='\t')
    df[['x','y', 'YEAR']].to_csv(xy_txt, sep='\t')

    return df


years = range(2002, 2013)
txt_template = '/Users/SamHooper/Grad_School/Thesis/data/eBird_data/erd_us48_data_grouped_by_year_v5.0/%s/checklists.csv'
states = ['Oregon', 'Washington', 'California']
species = 'Pheucticus_melanocephalus'
columns = {'SAMPLING_EVENT_ID': object, 'LATITUDE': float, 'LONGITUDE': float, 'YEAR': int, 'MONTH': int, 'DAY': int,
           'TIME': float, 'STATE_PROVINCE': object, 'COUNT_TYPE': object, 'EFFORT_HRS': float, 'EFFORT_D'
           'PRIMARY_CHECKLIST_FLAG': int}
state_field = 'STATE_PROVINCE'
nan_cols = ['LATITUDE', 'LONGITUDE', 'YEAR', 'DAY', 'TIME', 'COUNT_TYPE', 'EFFORT_HRS', species]
dup_col = 'PRIMARY_CHECKLIST_FLAG'
src_shp = '/Users/SamHooper/Grad_School/Thesis/data/eBird_data/erd_by_species/shp/WA_OR_CA_wgs84.shp'
dst_shp = '/Users/SamHooper/Grad_School/2015_2Spring/GEO584_spatial_stats/general_data/WA_OR_CA_albers.shp'
out_dir = '/Users/SamHooper/Grad_School/Thesis/data/eBird_data/erd_by_species'
xy_txt = '/Users/SamHooper/Grad_School/Thesis/data/eBird_data/erd_by_species/xy_prj.txt'
drop_cols = ['LATITUDE', 'LONGITUDE', state_field, dup_col, species, 'SAMPLING_EVENT_ID', 'MONTH']

# txt = os.path.join(out_dir, 'test.txt')
lat_field, lon_field = 'LATITUDE', 'LONGITUDE'

main(years, txt_template, states, species, columns, state_field, lat_field, lon_field, src_shp, dst_shp,
    out_dir, xy_txt=None, nan_cols=nan_cols, dup_col=dup_col, drop_cols=drop_cols)
