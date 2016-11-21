# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:20:49 2016

@author: shooper
"""

    if min_val == None or max_val == None:
        min_val = np.min(ar)
        max_val = np.max(ar)
    ds = None
    if min_val == nodata:
        min_val = min_val + 1
    # Calculate number of samples per bin and bin width
    #data_range = max_val - min_val
    samples_per = n_samples/len(bins)# + 1)
    #step = data_range/n_bins #increase per bin
    #bins = [(i * step, i * step + step) for i in range(n_bins)]
    #bins.append(nodata)
    
    print 'Making arrays of row and col indices... %s\n' % datetime.now()
    # Make 2D arrays of row and column indices
    ar_rows, ar_cols = np.indices(ar.shape)
    
    # Get random samples for each bin
    rows = []
    cols = []
    '''while i < n_bins:
    this_max = min_val + step
    print 'Getting random samples between %s and %s...' % (min_val, this_max)
    print datetime.now(), '\n'
    mask = (ar >= min_val) & (ar < this_max) & (ar != nodata)
    #these_rows, these_cols = np.where((ar >= min_val) & (ar < this_max) & (ar != nodata))'''
    nodata_mask = ar != nodata
    for b in bins:
        '''if b == nodata:
            print 'Getting random samples for nodata value %s...' % nodata
            print datetime.now(), '\n'
            mask = ar == nodata
            these_rows = ar_rows[mask]
            these_cols = ar_cols[mask]