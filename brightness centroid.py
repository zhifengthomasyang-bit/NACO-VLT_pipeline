# -*- coding: utf-8 -*-
"""
Created on Fri May 29 17:17:13 2026

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


bright_center = []

with fits.open(r"D:\headers\alignment\calib 2004-01-12 Imaging expt-0.5 OP3 FULL OP6 Ks HS DR S13\2004-01-13T03_09_13.954 - 2004-01-13T03_15_12.695\2004-01-13T03_09_13.954 - 2004-01-13T03_15_12.695.fits") as hdul:
    data = hdul[0].data.astype(np.float32)
    bright_center.append(data)
    
bright_center = np.array(bright_center)

xc, yc = 75, 75
R = 30
y, x = np.indices(data.shape)
r = np.sqrt((x - xc)**2 + (y - yc)**2)

mask = r <= R
img = data.copy()
img[~mask] = 0
img[img < 0] = 0

total_flux = np.sum(img)
x_centroid = np.sum(x * img) / total_flux
y_centroid = np.sum(y * img) / total_flux

print("Brightness centroid:")
print(f"x = {x_centroid:.2f}")
print(f"y = {y_centroid:.2f}")