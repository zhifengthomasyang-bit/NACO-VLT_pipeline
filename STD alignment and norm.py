# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 11:45:00 2026

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

alignment_dir = r"D:\headers\calib\2005-05-09\STD\BD20_2680\expt-0.2 OP3 FULL OP6 L_prime HWD Un L27"
alignment_files = glob.glob(os.path.join(alignment_dir, "*.fits"))
alignment = []

output_dir = r"D:\headers\cutouts"

name = alignment_dir.split("\\calib\\")[-1]
name = name.replace("\\", " ")
filename = f"PSF {name}.fits"

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
        
        tmp = img.copy()
        tmp[~valid] = -np.inf

        top_idx = np.argpartition(tmp.ravel(), -35)[-35:]

        top_idx = top_idx[np.argsort(tmp.ravel()[top_idx])[::-1]]

        rows, cols = np.unravel_index(top_idx, img.shape)
        
        y_med = np.median(rows)
        x_med = np.median(cols)
        
        x0,y0 = x_med,y_med
        half_size = 40// 2
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
            f"PSF_40x40_{basename}"
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

print(cutout.shape)


# median_stack = np.nanmedian(cutout, axis=0)

# psf = median_stack.copy()
# psf -= np.median(psf)

# psf_cleaned = np.maximum(psf, 0)

# psf_norm = psf_cleaned / np.sum(psf_cleaned)

psf = np.nanmean(cutout, axis=0)
ny, nx = psf.shape
edge = np.hstack([
    psf[:5, :].ravel(),
    psf[-5:, :].ravel(),
    psf[:, :5].ravel(),
    psf[:, -5:].ravel()
])
bg = np.median(edge)

psf = psf - bg

psf_norm = psf / np.sum(psf)

hdu = fits.PrimaryHDU(psf_norm)

output_path = os.path.join(cutout_dir, filename)

hdu.writeto(output_path, overwrite=True)

print("Saved:", output_path)

vals = psf_norm.flatten()
plt.figure(figsize=(6,4))
plt.hist(vals, bins=100, log=True)
plt.xlabel("PSF value (normalized intensity)")
plt.ylabel("Number of pixels")
plt.title("PSF Pixel Value Distribution")
plt.show()