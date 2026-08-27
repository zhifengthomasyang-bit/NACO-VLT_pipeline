# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:36:19 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from numpy.lib.stride_tricks import sliding_window_view

with fits.open(r'D:\headers\calib\2004-01-12\Imaging\epxt-0.4 OP3 FULL OP4 J HS DR S13\Calib NACO.2004-01-13T03_01_51.396.fits') as hdul:
    data = hdul[0].data

img = np.squeeze(data).astype(float)

y0, x0 = 500,600
half_size = 150// 2

y1, y2 = y0 - half_size, y0 + half_size
x1, x2 = x0 - half_size, x0 + half_size 

cutout = img[y1:y2, x1:x2]

mask = np.isfinite(cutout)

ny, nx = cutout.shape
y, x = np.mgrid[:ny, :nx]

amplitude_init = np.nanmax(cutout) - np.nanmedian(cutout)

x_mean_init = nx / 2
y_mean_init = ny / 2

sigma_init = 3.0

offset = np.nanmedian(cutout)

g_init = (
    models.Gaussian2D(
        amplitude=amplitude_init,
        x_mean=x_mean_init,
        y_mean=y_mean_init,
        x_stddev=sigma_init,
        y_stddev=sigma_init,
        theta=0
    )
    + models.Const2D(offset)
)

fitter = fitting.LevMarLSQFitter()

g_fit = fitter(
    g_init,
    x[mask],
    y[mask],
    cutout[mask]
)

gauss = g_fit[0]

x_global = x1 + gauss.x_mean.value
y_global = y1 + gauss.y_mean.value

print(f"\nGlobal position:")
print(f"x = {x_global:.3f}")
print(f"y = {y_global:.3f}")

pad = np.pad(cutout, 3, mode='constant', constant_values=np.nan)
windows = sliding_window_view(pad, (7, 7))

mean_map = np.nanmean(windows, axis=(2, 3))

cutout_filled = cutout.copy()
cutout_filled[np.isnan(cutout)] = mean_map[np.isnan(cutout)]

img_filled = cutout_filled

ny, nx = img_filled.shape
y, x = np.mgrid[:ny, :nx]

flux = np.sum(img_filled)

xc = np.sum(x * img_filled) / flux
yc = np.sum(y * img_filled) / flux

x_global = x1 + xc
y_global = y1 + yc

print(f"\nCentroid position:")
print(f"x = {x_global:.3f}")
print(f"y = {y_global:.3f}")

xc = int(round(x_global))
yc = int(round(y_global))

x1_new = xc - half_size
x2_new = xc + half_size 

y1_new = yc - half_size
y2_new = yc + half_size 


cutout_global = img[y1_new:y2_new, x1_new:x2_new]

print(cutout_global.shape)


hdu = fits.PrimaryHDU(cutout_global)

output_path = r'D:\headers\cutout_75x75 Calib NACO.2004-01-13T03_01_51.396.fits'

hdu.writeto(output_path, overwrite=True)

print(f"Saved to: {output_path}")



