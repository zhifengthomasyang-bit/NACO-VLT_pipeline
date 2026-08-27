# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 02:06:52 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import os
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt


path = r"D:\headers\calib\2004-01-12\Spec\expt-1.0 MASK Slit_86mas OP3 FULL OP4 Grism2 OP6 SK HS DR S54"
calib_files = glob.glob(os.path.join(path, "*.fits"))
calib = []

output_file = os.path.join(path,
    "2004-01-13T04_49_17.947 - 2004-01-13T04_59_06.387.fits"
)

x0 = 550
y0 = 630

all_spectra =[]

for f in calib_files:

    data = fits.getdata(f).astype(float)

    row = data[y0, :]

    valid = np.isfinite(row)

    row_valid = row[valid]
    x_valid = np.arange(len(row))[valid]

    max_idx = np.argsort(row_valid)[-5:]
    max_positions = x_valid[max_idx]

    min_idx = np.argsort(row_valid)[:5]
    min_positions = x_valid[min_idx]

    x2 = np.mean(max_positions)
    x1 = np.mean(min_positions)

    x1 = int(x1)
    x2 = int(x2)
    
    spec1 = data[:, x1]
    spec2 = data[:, x2]

    spectrum = (spec2 - spec1) / 2


    all_spectra.append(spectrum)

all_spectra = np.array(all_spectra)

mean_spectrum = np.nanmean(
    all_spectra,
    axis=0
)

hdu = fits.PrimaryHDU(mean_spectrum)
hdu.writeto(
    output_file,
    overwrite=True
)

y_pixel = np.arange(len(mean_spectrum))


plt.figure(figsize=(8,5))

plt.plot(
    y_pixel,
    mean_spectrum
)

plt.xlabel("Y pixel (dispersion direction)")
plt.ylabel("ADU")
plt.title("Extracted spectrum")

plt.grid(True)

plt.show()
