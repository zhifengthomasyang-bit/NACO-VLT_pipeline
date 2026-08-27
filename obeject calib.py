# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:42:27 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import os 

object_dir = r"D:\headers\catalogue\2005-05-09\STD\classify expt\expt-1.0\OP3 FULL OP6 Ks\S372-S_9137_VIS HD DR S27"
object_files = glob.glob(os.path.join(object_dir, "*.fits"))
object_calib = []

dark_file = [r"D:\headers\master\DARK\2005-05-09 MASK FLM_13\expt-1.0 HD DR S27.fits"]
master_dark = []

flat_file = [r"D:\headers\master\FLAT,LAMP\2005-05-09 MASK FLM_27\FLAT,LAMP expt-3.0 OP3 FULL OP6 Ks HD DR S27.fits"]
master_flat = []

flat_min = 0.85
flat_max = 1.15

for f in object_files:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    object_calib.append(data)
    
object_calib = np.array(object_calib)

for file in dark_file:
    with fits.open(file) as hdul:
        data = hdul[0].data.astype(np.float32)
        master_dark.append(data)
        
master_dark = np.array(master_dark)

for f in flat_file:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    data[(data < flat_min) | (data > flat_max)] = np.nan
    master_flat.append(data)

master_flat = np.array(master_flat)

correct_object = (object_calib - master_dark)/master_flat

input_name = os.path.basename(object_files[0])
#output_name = input_name.replace("%3A", "_")
#output_name = "Calib " + output_name


#hdu = fits.PrimaryHDU(correct_object)
#hdu.writeto(output_name, overwrite=True)

for i, f in enumerate(object_files):
    input_name = os.path.basename(f)
    output_name = "Calib " + input_name.replace("%3A", "_")

    hdu = fits.PrimaryHDU(correct_object[i])
    hdu.writeto(os.path.join(output_name), overwrite=True)
