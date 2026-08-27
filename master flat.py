# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 11:22:42 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob

fits_files = glob.glob(r"D:\headers\catalogue\2005-05-09\L_prime SKY/*.fits")
dark_file = [r"D:\headers\master\DARK\expt-0.2 HWD Un L27.fits"]

data_list = []
master_dark = []


for f in fits_files:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    data_list.append(data)
    
data_list = np.array(data_list)

for file in dark_file:
    with fits.open(file) as hdul:
        data = hdul[0].data.astype(np.float32)
        master_dark.append(data)
        
master_dark = np.array(master_dark)

flat_corrected = data_list - master_dark

flat_norm = []

for flat in flat_corrected:
    clipped = sigma_clip(flat, sigma=3).filled(np.nan)
    med = np.nanmedian(clipped)
    flat_norm.append(flat / med)

flat_norm = np.array(flat_norm)

outliers_clip = sigma_clip(flat_norm, sigma=3, axis=0).filled(np.nan)
master_flat = np.nanmedian(outliers_clip, axis=0)

hdu = fits.PrimaryHDU(master_flat)
hdu.writeto("master_flat.fits", overwrite=True)