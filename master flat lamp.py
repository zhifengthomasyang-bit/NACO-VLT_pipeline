# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:06:04 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob

fits_files = glob.glob(r"D:\headers\catalogue\2005-05-09\FLAT,LAMP\classify expt\expt-15.0\OP3 FULL OP6 NB_1.64 HD DR S13\FLAT.LAMP/*.fits")
dark_file = [r"D:\headers\master\DARK,LAMP\2005-05-09\expt-15.0 OP3 FULL OP6 NB_1.64 HD DR S13.fits"]

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
    med = np.nanmedian(flat)
    flat_norm.append(flat / med)

flat_norm = np.array(flat_norm)

master_flat = np.nanmedian(flat_norm, axis=0)

hdu = fits.PrimaryHDU(master_flat)
hdu.writeto("FLAT,LAMP expt-15.0 OP3 FULL OP6 NB_1.64 HD DR S13.fits", overwrite=True)