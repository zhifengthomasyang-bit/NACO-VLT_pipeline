# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 01:12:24 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import os

path = r"D:\headers\catalogue\2004-01-12\OBJECT\classify expt\expt-1.0 MASK Slit_86mas OP3 FULL OP4 Grism2 OP6 SK HS DR S54"
A_files = glob.glob(os.path.join(path, "A", "*.fits"))
B_files = glob.glob(os.path.join(path, "B", "*.fits"))
A_calib = []
B_calib = []


flat_file = [r"D:\headers\master\SPEC,FLAT\2004-01-12\SPEC FLAT expt-5.0 OP3 FULL OP4 Grism2 OP6 SK HS DR S54.fits"]
master_flat = []

output_dir = r"D:\headers\calib\2004-01-12\Spec\expt-1.0 MASK Slit_86mas OP3 FULL OP4 Grism2 OP6 SK HS DR S54"


flat_min = 0.7
flat_max = 1.3

for f in A_files:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    A_calib.append(data)
    
A_calib = np.array(A_calib)

for f in B_files:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    B_calib.append(data)
    
B_calib = np.array(B_calib)


for f in flat_file:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    data[(data < flat_min) | (data > flat_max)] = np.nan
    master_flat.append(data)

master_flat = np.array(master_flat)

master_flat = master_flat[0]

corrected_AB = []

for A, B in zip(A_calib, B_calib):

    diff = A - B
    corrected = diff / master_flat
    corrected_AB.append(corrected)

corrected_AB = np.array(corrected_AB)

for i, (A_file, B_file) in enumerate(zip(A_files, B_files)):

    A_name = os.path.basename(A_file)
    B_name = os.path.basename(B_file)

    output_name = (
        "Calib "
        + A_name.replace("%3A", "_")
        + " - "
        + B_name.replace("%3A", "_")
    )

    hdu = fits.PrimaryHDU(corrected_AB[i])

    hdu.writeto(
        os.path.join(output_dir, output_name),
        overwrite=True
    )
