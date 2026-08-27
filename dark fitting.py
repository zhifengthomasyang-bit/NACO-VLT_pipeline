# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 13:02:37 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np

dark_files = [r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-0.4 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-0.5 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-0.7 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-1.0 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-2.0 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-3.5 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-4.0 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-9.0 HS DR L27.fits',
              r'D:\headers\master\DARK\2004-01-12 MASK FLM_27\expt-30.0 HS DR L27.fits']
exposure_times = np.array([0.4, 0.5, 0.7, 1.0, 2.0, 3.5, 4.0, 9.0, 30.0], dtype=np.float32)

def generate_pixel_models(files, times):
    master_dark = []
    for file in files:
        with fits.open(file) as hdul:
            data = hdul[0].data.astype(np.float32)
            master_dark.append(data)
    
    master_dark = np.array(master_dark)
    num_exposures, height, width = master_dark.shape
    
    A = np.vstack([times, np.ones(len(times))]).T
    flat_dark = master_dark.reshape(num_exposures, -1)
    params, residuals, rank, s = np.linalg.lstsq(A, flat_dark, rcond=None)
    k_map = params[0].reshape(height, width)
    b_map = params[1].reshape(height, width)
    res_map = np.sqrt(residuals / num_exposures).reshape(height, width)
    return k_map, b_map, res_map

k_map, b_map, res_map = generate_pixel_models(dark_files, exposure_times)
fits.writeto('pixel_slope_k.fits', k_map, overwrite=True)
fits.writeto('pixel_intercept_b.fits', b_map, overwrite=True)
fits.writeto('pixel_residuals.fits', res_map, overwrite=True)

def predict_dark(t, k, b):
    return k * t + b

exptime = 3.0
predicted_15s = predict_dark(exptime, k_map, b_map)
fits.writeto('expt-3.0 HS DR.fits', predicted_15s, overwrite=True)
