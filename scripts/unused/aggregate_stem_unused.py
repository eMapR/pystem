# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 16:18:30 2016

@author: shooper
"""

"""def fill_tile_band(tile_size, tile_coords, set_coords, ar_pred, tx, nodata):
    '''
    Fill an array of zeros of shape tile_size, located at tile_coords with an 
    offset array, ar_pred, located at set_coords
    '''
    # Calc offsets
    row_off, col_off = calc_offset(tile_coords[['ul_x','ul_y']],
                                   set_coords[['ul_x','ul_y']], tx)
    # Get the offset indices of each array 
    tile_inds, set_inds = mosaic.get_offset_array_indices(
        tile_size, ar_pred.shape, (row_off, col_off))
    tile_row_u, tile_row_d, tile_col_l, tile_col_r = tile_inds
    set_row_u,  set_row_d,  set_col_l,  set_col_r  = set_inds
    
    # Fill just the part of the array that overlaps
    ar_tile = np.full(tile_size, np.nan)
    ar_pred = ar_pred.astype(float)
    ar_pred[ar_pred == nodata] = np.nan
    try:
        ar_tile[tile_row_u:tile_row_d, tile_col_l:tile_col_r] =\
        ar_pred[set_row_u:set_row_d, set_col_l:set_col_r]
    except Exception as e:
        print e
        print '\nProblem with offsets'
        print row_off, col_off, set_coords, tile_coords      

    return ar_tile"""
    
def fill_tile_band(tile_size, tile_coords, set_coords, ar_pred, tx, nodata):
    '''
    Fill an array of zeros of shape tile_size, located at tile_coords with an 
    offset array, ar_pred, located at set_coords
    '''
    # Calc offsets
    row_off, col_off = calc_offset_from_tile(tile_coords[['ul_x','ul_y']],
                                             set_coords[['ul_x','ul_y']],
                                             tx)
    # Get the offset indices of both arrays 
    tile_inds, set_inds = mosaic.get_offset_array_indices(
        tile_size, ar_pred.shape, (row_off, col_off))
    tile_row_u, tile_row_d, tile_col_l, tile_col_r = tile_inds
    set_row_u,  set_row_d,  set_col_l,  set_col_r  = set_inds
    
    # Fill just the part of the array that overlaps
    ar_tile = np.full(tile_size, np.nan)
    ar_pred = ar_pred.astype(float)
    ar_pred[ar_pred == nodata] = np.nan
    try:
        ar_tile[tile_row_u:tile_row_d, tile_col_l:tile_col_r] =\
        ar_pred[set_row_u:set_row_d, set_col_l:set_col_r]
    except Exception as e:
        print e
        print '\nProblem with offsets'
        print row_off, col_off, set_coords, tile_coords      

    return ar_tile