# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 21:06:28 2016

@author: shooper
"""

def plot_sets_on_shp(ds_coords, ds_extent, max_size, df_sets, support_size, out_dir=None, fill='0.7', line_color='0.5', line_width=1.0, pad=.54):
    '''
    Return an RGB array of a shapefile. 
    '''
    
    # Get dataset size and compute ratio to rows and cols
    x_min, x_max, y_min, y_max = ds_extent
    delta_x = x_max - x_min
    delta_y = y_max - y_min
    
    # Figure out which dimension is larger, make it max_size, and calculate
    #    the proportional size of the other dimension
    if delta_x >= delta_y: 
        cols = max_size
        rows = int(max_size * delta_y/delta_x)
    else:
        cols = int(max_size * delta_x/delta_y)
        rows = max_size
    
    # Compute scale factors to calculate rows and columns from coords
    x_scale = 1.0/delta_x
    y_scale = 1.0/delta_y
    
    # Create the plot
    fig = plt.figure(figsize=(cols/72.0, rows/72.0), dpi=72)#, tight_layout={'pad':pad})
    sub = fig.add_subplot(1, 1, 1, axisbg='w', frame_on=False)
    sub.axes.get_yaxis().set_visible(False)
    sub.axes.get_xaxis().set_visible(False) 
    #sub.set_xlim([x_min, x_max])
    #sub.set_ylim([y_min, y_max])
    
    # Make a list of patches where each feature is a separate patch 
    patches = []
    ds_coords = [c for c in ds_coords if c != None]
    for feature in ds_coords:
        img_coords = [((pt[0] - x_min) * x_scale, (pt[1] - y_min) * y_scale) for pt in feature]
        print img_coords
        print ''
        #img_coords = [(pt[0], pt[1]) for pt in feature]
        ar_coords = np.array(img_coords)
        poly = matplotlib.patches.Polygon(ar_coords)
        patches.append(poly)
        #sub.add_patch(poly, facecolor=fill, lw=line_width, edgecolor=line_color)
    
    # Make the patch collection and add it to the plot
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, color=fill, lw=line_width, edgecolor=line_color)
    sub.add_collection(p)
    
    # Plot the support sets
    for ind, r in df_sets.iterrows():
        lr_x = (r['min_x'] - x_min) * x_scale
        lr_y = (r['min_y'] - y_min) * y_scale
        w = support_size[1] * x_scale
        h = support_size[0] * y_scale
        #if lr_x <.2 and lr_y<.2:
        sub.add_patch(plt.Rectangle((lr_x, lr_y), w, h, facecolor='none', lw='0.3', alpha=.2))
    plt.show()

def sample_support_set(df, pct_train):
    '''
    Return train and test dfs from random samples such that train and test
    don't contain any of the same locations
    '''
    n_train = int(len(df) * pct_train)
    
    # Randomly sample undf up to n_training samples with replacement
    try:
        inds = rnd.sample(df.index, n_train)
        df_train = df.ix[inds]
        df_test = df.drop(inds) 
    except ValueError:
        # A value error likely means n_train > len(df)
        return None, None
    
    # Drop any training points that share the same location. ~ works as !
    df_train = df_train[~df_train['row'].isin(df_test['row']) &
    ~df_train['col'].isin(df_test['col'])]
    
    return df_train, df_test
    

#lr = [[x + (x_size * x_res_sign), y + (y_size * y_res_sign)] for x, y in ul]
#df = pd.DataFrame(cells, columns=['ul_x', 'ul_y', 'lr_x', 'lr_y'])

def sample_support_set(df, pct_train):
    '''
    Return train and test dfs from random samples such that train and test
    don't contain any of the same locations
    '''
    # Get a list of unique row,col pairs
    unique = [(rc[0],rc[1]) for rc in df[['row','col']].drop_duplicates().values]
    n_train = int(len(unique) * pct_train)
    
    # Randomly sample n_training locations without replacement
    try:
        rowcol = random.sample(unique, n_train) #Sample unique row,col
        rows = [rc[0] for rc in rowcol]
        cols = [rc[1] for rc in rowcol]
        # Get all of the obs from those row cols for training points
        df_train = df[df.row.isin(rows) & df.col.isin(cols)]
        # Get all of the obs not from those row,cols for testing
        #df_test = df[~df['row'].isin(rows) & ~df['col'].isin(cols)]

    except ValueError:
        # A value error likely means n_train > len(df). This isn't
        #   really possible with the current code, but whatever.
        return None, None
    
    return df_train#, df_test  