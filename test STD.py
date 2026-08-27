# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 11:31:34 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from numpy.lib.stride_tricks import sliding_window_view


with fits.open(r'D:\headers\calib\2004-01-12\STD\S165-E_9136_VIS\expt-0.5 OP3 FULL OP4 J HS DR S27\Calib NACO.2004-01-13T06_44_16.326.fits') as hdul:
    data = hdul[0].data

img = np.squeeze(data).astype(float)

valid = np.isfinite(img)

values = img[valid]

tmp = img.copy()
tmp[~valid] = -np.inf

top_idx = np.argpartition(tmp.ravel(), -10)[-10:]

top_idx = top_idx[np.argsort(tmp.ravel()[top_idx])[::-1]]

rows, cols = np.unravel_index(top_idx, img.shape)

print("Top 10 brightest pixels:")
for i, (r, c) in enumerate(zip(rows, cols), 1):
    print(
        f"{i}: value={img[r,c]:.2f}, "
        f"(y={r}, x={c})"
    )
    
y_med = np.median(rows)
x_med = np.median(cols)
print(f"\nMedian coordinate:")
print(f"(y={y_med:.1f}, x={x_med:.1f})")

x0,y0 = x_med,y_med
half_size = 20// 2
x0 = int(x0)
y0 = int(y0)

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

hdu = fits.PrimaryHDU(cutout)

out_path = r"D:\psf_cutout_20x20.fits"
hdu.writeto(out_path, overwrite=True)

print("Saved:", out_path)