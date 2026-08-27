# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:44:15 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import os

path = r"D:\headers\catalogue\2004-01-12\WAVE,LAMP\SK"
on_files = glob.glob(os.path.join(path, "ON", "*.fits"))
off_files = glob.glob(os.path.join(path, "OFF", "*.fits"))

file_name = path.split('WAVE,LAMP\\')[-1].replace('\\', ' ') + ".fits"


dark_stack = []

for file in off_files:
    with fits.open(file) as hdul:
        data = hdul[0].data.astype(np.float32)
        dark_stack.append(data)

dark_stack = np.array(dark_stack)


master_dark = np.median(dark_stack, axis=0)

output_path = os.path.join(path, file_name)

hdu = fits.PrimaryHDU(master_dark)
hdu.writeto(output_path, overwrite=True)


data_list = []
for f in on_files:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    data_list.append(data)
    
data_list = np.array(data_list)

    
flat_corrected = data_list - master_dark

master_flat_corrected = np.median(flat_corrected, axis=0)
output_master_path = os.path.join(path, f"FLAT {file_name}")

hdu = fits.PrimaryHDU(master_flat_corrected.astype(np.float32))
hdu.writeto(output_master_path, overwrite=True)