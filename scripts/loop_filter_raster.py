# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 22:08:02 2016

@author: shooper
"""
import os
import sys
import filter_raster

def main(param_dir):
    param_files = [os.path.join(param_dir, f) for f in os.listdir(param_dir)]
    n_files = len(param_files)
    for i, params in enumerate(param_files):
        print 'Filtering %s of %s rasters: %s' % (i + 1, n_files, os.path.basename(params))
        filter_raster.main(params)
        

if __name__ == '__main__':
    param_dir = sys.argv[1]
    sys.exit(main(param_dir))