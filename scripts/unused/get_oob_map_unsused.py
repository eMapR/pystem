# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:28:28 2016

@author: shooper
"""

                '''test_samples = this_oob[val_col]
                test_predictors = this_oob[var_cols]

                # Get oob score and fill tile with that val
                #try:
                with open(s_row.dt_file, 'rb') as f: 
                    dt_model = pickle.load(f)
                ar_oob_rate, tile_inds, oob_rate = get_oob_array(dt_model, test_samples, test_predictors, tx, tile_ul, tile_size, set_coords, err_threshold)#'''
                
                
def get_oob_array(dt, test_samples, test_predictors, tx, tile_ul, tile_size, set_coords, err_threshold):
    
    ''' rewrite so that I calculate all oob rates and store in df sets, then create the arrays on the fly'''
    # get OOB score from dt
    oob_rate = stem.calc_oob_rate(dt, test_samples, test_predictors, err_threshold)
    
    # Calc array size
    set_ulx, set_uly, set_lrx, set_lry = set_coords
    offset = stem.calc_offset(tile_ul, (set_ulx, set_uly), tx)
    xsize = abs((set_ulx - set_lrx)/tx[1])
    ysize = abs((set_uly - set_lry)/tx[5])
    t_inds, a_inds = mosaic.get_offset_array_indices(tile_size, (ysize, xsize), offset)
    nrows = a_inds[1] - a_inds[0]
    ncols = a_inds[3] - a_inds[2]
    
    ar = np.full((nrows, ncols), oob_rate, dtype=np.int16)
    
    return ar, t_inds, oob_rate