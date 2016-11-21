# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 22:57:24 2016

@author: shooper
"""
import os
import sys
import fnmatch
import predict_stem

def main(lower_yr, upper_yr, target, param_dir):
    years = range(lower_yr, upper_yr + 1)
    files = os.listdir(param_dir)
    for year in years[: : -1]:
        print 'Predicting for %s for year %s' % (target, year)
        try:
            params_fn = fnmatch.filter(files, '*params*%s%s.txt' % (target, year))[0]
        except:
            import pdb; pdb.set_trace()
            continue
        predict_stem.main(os.path.join(param_dir, params_fn))

if __name__ == '__main__':
    lower_yr, upper_yr, target, param_dir = sys.argv[1:]
    lower_yr = int(lower_yr)
    upper_yr = int(upper_yr)
    sys.exit(main(lower_yr, upper_yr, target, param_dir))