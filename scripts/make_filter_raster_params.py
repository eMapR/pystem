# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 20:53:21 2016

@author: shooper
"""

import os


def write_area_params(raster_path, years, res, vals, n_jobs, out_dir):
    for val in vals:
        name = vals[val]
        for r in res:
            stamp = '%s_%s' % (name, r)
            for year in years:
                out_bn = 'filter_raster_params_%s_%s.txt' % (stamp, year)
                out_txt = os.path.join(out_dir, out_bn)
                path = raster_path % year
                out_path = '/vol/v2/stem/ebird/predictors/{0}/{0}_{1}.bsq'.format(stamp, year)
                print name, r, year
                with open(out_txt, 'w') as txt:
                    txt.write('Parameters for filter_raster.py\n')
                    txt.write('\n')
                    txt.write('path; %s\n' % path)
                    txt.write('databand; 1\n')
                    txt.write('function; area\n')
                    txt.write('kernel_size; %s\n' % (int(r.replace('m', ''))/15))
                    txt.write('out_path; %s\n' % out_path)
                    txt.write('\n')
                    txt.write('Optional parameters\n')
                    txt.write('n_tiles; 25, 15\n')
                    txt.write('n_jobs; %s\n' % n_jobs)
                    txt.write('nodata; 0\n')
                    txt.write('kernel_type; circle\n')
                    txt.write('out_nodata; 255\n')
                    txt.write('filter_value; %s\n' % val)


def write_avg_params(raster_path, years, res, names, n_jobs, out_dir):

    for name in names:
        for r in res:
            stamp = 'avg%s_%s' % (name, r)
            for year in years:
                out_bn = 'filter_raster_params_%s_%s.txt' % (stamp, year)
                out_txt = os.path.join(out_dir, out_bn)
                path = raster_path.format(name, year)
                out_path = '/vol/v2/stem/ebird/predictors/{0}/{0}_{1}.bsq'.format(stamp, year)
                print name, r, year
                with open(out_txt, 'w') as txt:
                    txt.write('Parameters for filter_raster.py\n')
                    txt.write('\n')
                    txt.write('path; %s\n' % path)
                    txt.write('databand; 1\n')
                    txt.write('function; average\n')
                    txt.write('kernel_size; %s\n' % (int(r.replace('m', ''))/15))
                    txt.write('out_path; %s\n' % out_path)
                    txt.write('\n')
                    txt.write('Optional parameters\n')
                    txt.write('n_tiles; 25, 15\n')
                    txt.write('n_jobs; %s\n' % n_jobs)
                    txt.write('nodata; 255\n')
                    txt.write('kernel_type; circle\n')




out_dir = '/vol/v2/stem/ebird/predictors/param_files'
years = range(2002, 2008)[::-1]
res = ['225m']
n_jobs = 30
vals = {11: 'pctwater', 12: 'pctsnow', 23: 'pcturban', 31: 'pctbareground', 41: 'pctdeciduous', 42: 'pctconiferous', 52: 'pctshrub'}
raster_path = '/vol/v1/proj/lst/outputs/models/randomforest/rfprediction_mosaic/yearly/lst_run1_prediction_voting_lulc_RF_mosaic_%s.bsq'
#write_area_params(raster_path, years, res, vals, n_jobs, out_dir)

raster_path = '/vol/v2/stem/{0}/time_series/{1}/{0}_vote_{1}.bsq'   
names = ['imperv']  
write_avg_params(raster_path, years, res, names, n_jobs, out_dir)         
                
    