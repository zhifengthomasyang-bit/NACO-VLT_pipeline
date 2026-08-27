# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:20:31 2026

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

alignment_dir = r"D:\headers\calib\2005-05-09\Imaging\Pallas9\expt-2.0 OP3 FULL OP6 Ks HD DR S13"
alignment_files = glob.glob(os.path.join(alignment_dir, "*.fits"))
alignment = []

output_dir = r"D:\headers\cutouts"

os.makedirs(output_dir, exist_ok=True)



for f in alignment_files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(np.float32)
        alignment.append(data)
        
alignment = np.array(alignment)

fitter = fitting.LevMarLSQFitter()

for f in alignment_files:

    print(f"\nProcessing: {f}")
    with fits.open(f) as hdul:
        data = hdul[0].data

    img = np.squeeze(data).astype(float)
      
    valid = np.isfinite(img)

    values = img[valid]
    
    threshold = np.percentile(values, 99.9)

    mask = img >= threshold

    y, x = np.indices(img.shape)

    x_sel = x[mask]
    y_sel = y[mask]
    w = img[mask]

    x_centroid = np.sum(x_sel * w) / np.sum(w)
    y_centroid = np.sum(y_sel * w) / np.sum(w)
    
    x0 = int(x_centroid)
    y0 = int(y_centroid)
    
    row_profile = img[y0, :]
    col_profile = img[:, x0]
    
    # plt.subplot(1,2,1)
    # plt.plot(row_profile)
    # plt.title(f"Row Profile (y={y0})")
    # plt.xlabel("X pixel")
    # plt.ylabel("Intensity")


    # plt.subplot(1,2,2)
    # plt.plot(col_profile)
    # plt.title(f"Column Profile (x={x0})")
    # plt.xlabel("Y pixel")
    # plt.ylabel("Intensity")

    # plt.tight_layout()
    # plt.show()

    x_row = np.arange(len(row_profile))
    x_col = np.arange(len(col_profile))

    mask_row = np.isfinite(row_profile) & (row_profile > 0)
    mask_col = np.isfinite(col_profile) & (col_profile > 0)

    x_row_fit = x_row[mask_row]
    y_row_fit = row_profile[mask_row]

    x_col_fit = x_col[mask_col]
    y_col_fit = col_profile[mask_col]

    background_row = np.median(y_row_fit)
    background_col = np.median(y_col_fit)
    
    g_init_row = (models.Gaussian1D(
        amplitude=np.max(y_row_fit),
        mean=x_row_fit[np.argmax(y_row_fit)],
        stddev=5
    )+
    models.Const1D(background_row)
    )

    g_fit_row = fitter(
        g_init_row,
        x_row_fit,
        y_row_fit
    )

    g_init_col = (models.Gaussian1D(
        amplitude=np.max(y_col_fit),
        mean=x_col_fit[np.argmax(y_col_fit)],
        stddev=5
    )+models.Const1D(background_col)
    )

    g_fit_col = fitter(
        g_init_col,
        x_col_fit,
        y_col_fit
    )

    row_center = g_fit_row[0].mean.value
    col_center = g_fit_col[0].mean.value


    x0 = int(round(row_center))
    y0 = int(round(col_center))
    


    half_size = 150 // 2

    y1 = y0 - half_size
    y2 = y0 + half_size

    x1 = x0 - half_size
    x2 = x0 + half_size

    cutout = img[y1:y2, x1:x2]


    mask = np.isfinite(cutout)

    ny, nx = cutout.shape

    y, x = np.mgrid[:ny, :nx]

    amplitude_init = (
        np.nanmax(cutout)
        - np.nanmedian(cutout)
    )

    offset = np.nanmedian(cutout)

    g_init = (
        models.Gaussian2D(
            amplitude=amplitude_init,
            x_mean=nx / 2,
            y_mean=ny / 2,
            x_stddev=3,
            y_stddev=3,
            theta=0
        )
        + models.Const2D(offset)
    )

    g_fit = fitter(
        g_init,
        x[mask],
        y[mask],
        cutout[mask]
    )

    gauss = g_fit[0]

    x_global = x1 + gauss.x_mean.value
    y_global = y1 + gauss.y_mean.value


    xc = int(round(x_global))
    yc = int(round(y_global))

    x1_new = xc - half_size
    x2_new = xc + half_size

    y1_new = yc - half_size
    y2_new = yc + half_size

    cutout_global = img[
        y1_new:y2_new,
        x1_new:x2_new
    ]

    basename = os.path.basename(f)

    output_name = (
        f"cutout_150x150_{basename}"
    )

    output_path = os.path.join(
        output_dir,
        output_name
    )


    hdu = fits.PrimaryHDU(cutout_global)

    hdu.writeto(
        output_path,
        overwrite=True
    )

    print("Saved:", output_path)

cutout_dir = r"D:\headers\cutouts"
cutout_files = glob.glob(os.path.join(cutout_dir, "*.fits"))
cutout = []

for f in cutout_files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(np.float32)
        cutout.append(data)
    
cutout = np.array(cutout)

median_stack = np.nanmedian(cutout, axis=0)

hdu = fits.PrimaryHDU(median_stack)

output_path = os.path.join(cutout_dir, "calib 2005-05-09 Imaging Pallas9 expt-2.0 OP3 FULL OP6 Ks HD DR S13.fits")

hdu.writeto(output_path, overwrite=True)

print("Saved:", output_path)

