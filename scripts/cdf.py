# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 21:44:25 2016

@author: shooper
"""
import gdal
import numpy as np
from matplotlib import pyplot as plt

plt.style.use('ggplot')

def cdf(ar, bins):
    
    c, b = np.histogram(ar, bins=bins)
    cdf = [c[:i + 1].sum() for i in b[:-1]]
    
    return cdf


def main(p_path, t_path, p_mean, bins, xlim=None, ylim=None, out_png=None):
    
    ds = gdal.Open(p_path)
    ar_p = ds.ReadAsArray()
    cdf_p = cdf(ar_p, bins)
    ds = None
    
    ds = gdal.Open(t_path)
    ar_t = ds.ReadAsArray()
    cdf_t = cdf(ar_t, bins)
    ds = None
    
    ds = gdal.Open(p_mean)
    ar_m = ds.ReadAsArray()
    cdf_m = cdf(ar_m, bins)
    
    bins = bins[:-1]
    plt.plot(bins, cdf_p, '-', label='STEM mode')
    plt.plot(bins, cdf_m, '-', label='STEM mean')
    plt.plot(bins, cdf_t, 'k--', label='NLCD')
    
    if xlim:
        plt.xlim(xlim[0], xlim[1])
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    plt.legend(loc='lower right', frameon=False)
    plt.tick_params(axis='x', width=0)
    plt.tick_params(axis='y', width=0)
    
    plt.xlabel('Map Value')
    plt.ylabel('Pixels (x 10 )')
    if out_png:
        plt.savefig(out_png)
    del ar_t, ar_p
    plt.clf()

p_vote = '/vol/v2/stem/canopy/time_series/2011/canopy_vote_2011.bsq'
p_mean = '/vol/v2/stem/canopy/time_series/2011/canopy_2011_mean.bsq'
t_path = '/vol/v2/stem/canopy/truth_map/canopy2011_CAORWA.tif'
bins = range(102)
xlim = [-10, 110]
ylim = [0, 10**9]
out_png = '/vol/v2/stem/canopy/time_series/2011/cdf_plot_mean_vote.png'
main(p_vote, t_path, p_mean, bins, xlim, None, out_png)

'''p_path = '/vol/v2/stem/imperv/time_series/2011/imperv_vote_2011.bsq'
t_path = '/vol/v2/stem/imperv/truth_map/imperv2011_CAORWA.tif'
out_png = '/vol/v2/stem/imperv/time_series/2011/cdf_plot.png'
main(p_path, t_path, bins, xlim, ylim, out_png)'''