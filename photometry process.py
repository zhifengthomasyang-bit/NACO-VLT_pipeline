# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 15:51:41 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from numpy.lib.stride_tricks import sliding_window_view
import os 

alignment_dir = r"D:\headers\alignment\calib 2005-05-09 Imaging Ceres9_L expt-2.0 OP3 FULL OP6 Ks HD DR S13\2"
alignment_files = glob.glob(os.path.join(alignment_dir, "*.fits"))
alignment = []

output_filename = "2005-05-10T01_18_33.644 - 2005-05-10T01_23_16.090.fits"
output_path = os.path.join(alignment_dir, output_filename)

expt = 2.0

for f in alignment_files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(np.float32)
        alignment.append(data)
        
alignment = np.array(alignment)

median_image = np.nanmedian(alignment, axis=0)

corner_size = 5

h, w = median_image.shape
corners = [
    median_image[0:corner_size, 0:corner_size],         
    median_image[0:corner_size, w-corner_size:w],        
    median_image[h-corner_size:h, 0:corner_size],       
    median_image[h-corner_size:h, w-corner_size:w]      
]

background = np.min(corners)

final_image = (median_image - background) / expt

hdu = fits.PrimaryHDU(final_image)
hdu.writeto(output_path, overwrite=True)