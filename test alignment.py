# -*- coding: utf-8 -*-
"""
Created on Tue May 26 10:00:12 2026

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

valid = np.isfinite(img)

values = img[valid]

threshold = np.percentile(values, 99.975)

mask = img >= threshold

y, x = np.indices(img.shape)

x_sel = x[mask]
y_sel = y[mask]
w = img[mask]

x_centroid = np.sum(x_sel * w) / np.sum(w)
y_centroid = np.sum(y_sel * w) / np.sum(w)

print("Centroid (top 0.3%)")
print("x =", x_centroid)
print("y =", y_centroid)

x0 = int(x_centroid)
y0 = int(y_centroid)

row_profile = img[y0, :]
col_profile = img[:, x0]
x_row = np.arange(len(row_profile))
x_col = np.arange(len(col_profile))

mask_row = np.isfinite(row_profile) & (row_profile > 0)

x_row_fit = x_row[mask_row]
y_row_fit = row_profile[mask_row]

mask_col = np.isfinite(col_profile) & (col_profile > 0)

x_col_fit = x_col[mask_col]
y_col_fit = col_profile[mask_col]

background_row = np.median(y_row_fit)
background_col = np.median(y_col_fit)

g_init_row = (models.Gaussian1D(
    amplitude=np.max(y_row_fit),
    mean=x_row_fit[np.argmax(y_row_fit)],
    stddev=5
)+models.Const1D(background_row)
)

g_init_col = (models.Gaussian1D(
    amplitude=np.max(y_col_fit),
    mean=x_col_fit[np.argmax(y_col_fit)],
    stddev=5
)+models.Const1D(background_col)
)

fitter = fitting.LevMarLSQFitter()

g_fit_row = fitter(
    g_init_row,
    x_row_fit,
    y_row_fit
)

g_fit_col = fitter(
    g_init_col,
    x_col_fit,
    y_col_fit
)

row_center = g_fit_row[0].mean.value
col_center = g_fit_col[0].mean.value

print("Row Gaussian center =", row_center)
print("Column Gaussian center =", col_center)


plt.subplot(1,2,1)
plt.plot(row_profile)
plt.title(f"Row Profile (y={y0})")
plt.xlabel("X pixel")
plt.ylabel("Intensity")


plt.subplot(1,2,2)
plt.plot(col_profile)
plt.title(f"Column Profile (x={x0})")
plt.xlabel("Y pixel")
plt.ylabel("Intensity")

plt.tight_layout()
plt.show()

x0,y0 = row_center,col_center
half_size = 150// 2
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

fig, axes = plt.subplots(1, 2, figsize=(12, 4))


axes[0].plot(x_row_fit, y_row_fit, label='Row data')

axes[0].plot(
    x_row_fit,
    g_fit_row(x_row_fit),
    label='Gaussian fit'
)

axes[0].axvline(
    row_center,
    color='r',
    linestyle='--',
    label=f'center = {row_center:.2f}'
)

axes[0].set_title('Row profile')
axes[0].set_xlabel('X pixel')
axes[0].set_ylabel('Flux')
axes[0].legend()


axes[1].plot(x_col_fit, y_col_fit, label='Column data')

axes[1].plot(
    x_col_fit,
    g_fit_col(x_col_fit),
    label='Gaussian fit'
)

axes[1].axvline(
    col_center,
    color='r',
    linestyle='--',
    label=f'center = {col_center:.2f}'
)

axes[1].set_title('Column profile')
axes[1].set_xlabel('Y pixel')
axes[1].set_ylabel('Flux')
axes[1].legend()

plt.tight_layout()
plt.show()

