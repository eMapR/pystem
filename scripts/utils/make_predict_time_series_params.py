#!/usr/bin/env python
"""
Make parameter text files for to create a time series from predict_stem.py.

Usage:
    make_predict_time_series_params.py <predict_params> <start_year> <end_year> <out_dir> <txt_out_dir> [--n_jobs=<int>] [--agg_stats=<str>] [--confusion=<b>] [--subset_shp=<path>] [--n_tiles=<str>]
    make_predict_time_series_params.py -h | --help

Options:
    -h --help               Show this screen.
    --n_jobs=<int>          Integer number of cores to use
    --agg_stats=<str>       list of agg stats[default: mean, vote, median, stdv]
    --confusion=<b>         Specifies whether to evaluate with a confusion matrix
    --subset_shp=<path>     If specified, AOI subset to process prediction within
    --n_tiles=<str>         'y, x' str of tiles for predicting
"""


import sys
import os
import re
import docopt
package_dir = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(package_dir)
from stem import read_params


def main(predict_params, start_year, end_year, out_dir, txt_out_dir, n_jobs=0, agg_stats='mean, vote, median, stdv', confusion=False, subset_shp=None, n_tiles=None):
    
    param_dict, df_var_orig = read_params(predict_params)
    for k, v in param_dict.iteritems():
        param_dict[k] = v.replace('"','')
    
    param_basename = os.path.basename(predict_params)
    out_txt = os.path.join(txt_out_dir, param_basename)
    if not os.path.isdir(txt_out_dir):
        os.mkdir(txt_out_dir)
    
    
    for year in range(int(start_year), int(end_year) + 1):
        print 'Making params for year %s...' % year
        
        df_var = df_var_orig.copy()     
        
        # Write the variable param table first
        band = year - 1983 # 1983 becuse gdal bands start at 1
        df_var.ix[df_var.data_band != 1, 'data_band'] = band
        df_var.data_band = df_var.data_band.astype(int)
        
        this_txt = out_txt.replace('.txt', '_%s.txt' % year)
        df_var.to_csv(this_txt, sep='\t')
        
        # Adjust a couple of variable values
        file_stamp = os.path.basename(param_dict['model_dir'].split('_')[0]) + '_' + str(year)
        param_dict['file_stamp'] = file_stamp
        param_dict['out_dir'] = os.path.abspath(os.path.join(out_dir, str(year)))
        param_dict['agg_stats'] = agg_stats
        #if 'confusion_params' in param_dict and not confusion:
            #del param_dict['confusion_params']
        
        with open(this_txt, 'a') as txt:
            txt.write('\n')
            txt.write('model_dir; %s\n' % param_dict['model_dir'])
            txt.write('train_params; %s\n' % param_dict['train_params'])
            txt.write('mosaic_path; %s\n' % param_dict['mosaic_path'])
            txt.write('support_size; %s\n' % param_dict['support_size'])
            txt.write('nodata; %s\n' % param_dict['nodata'])
            txt.write('out_dir; %s\n' % param_dict['out_dir'])
            txt.write('agg_stats; %s\n' % param_dict['agg_stats'])
            txt.write('\nOptional Parameters\n')
            txt.write('file_stamp; %s\n' % param_dict['file_stamp'])
            if int(n_jobs) != 0:
                n_jobs = int(n_jobs)
                txt.write('n_jobs; %s\n' % n_jobs)
            if subset_shp:
                txt.write('subset_shp; %s\n' % subset_shp)
            if n_tiles:
                txt.write('n_tiles; %s\n' % n_tiles)
        
        print 'Params written to %s\n' % this_txt

'''if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))#'''

if __name__ == '__main__':
    cl_args = docopt.docopt(__doc__)
    
    #except docopt.DocoptExit as e:
    #    print e.message
    
    # get rid of extra characters from doc string and 'help' entry
    args = {re.sub('[<>-]*', '', k): v for k, v in cl_args.iteritems()
            if k != '--help' and k != '-h'}
    
    sys.exit(main(**args))#'''
    
'''import glob
predict_params = '/vol/v2/stem/conus/models/imperv_20171026_1506/predict_stem_params.txt'
start_year = 1990
end_year = 2016
out_dir = '/vol/v2/stem/conus/time_series/imperv'
txt_out_dir = '/vol/v2/stem/conus/time_series/imperv/param_files'
names = 'detroit', 'sacramento', 'denver', 'dc', 'vegas', 'orlando'


for n in names:#glob.glob('/vol/v1/general_files/user_files/samh/pecora/imperv_maps/urban_extents/*.shp'):
    shp = '/vol/v1/general_files/user_files/samh/pecora/imperv_maps/urban_extents/%s.shp' % n
    main(predict_params, start_year, end_year, out_dir, txt_out_dir, 20, 0, subset_shp=shp)'''