# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 09:28:21 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

file_path = r"D:\headers\PSF\PSF 2005-05-09 STD BD20_2680 expt-0.2 OP3 FULL OP6 L_prime HWD Un L27\PSF 2005-05-09 STD BD20_2680 expt-0.2 OP3 FULL OP6 L_prime HWD Un L27.fits"
output = r"D:\headers\PSF\PSF 2005-05-09 STD BD20_2680 expt-0.2 OP3 FULL OP6 L_prime HWD Un L27 2005-05-10T08_47_01.833 - 08_52_13.278.png"

hdul = fits.open(file_path)

data = hdul[0].data

if data is None:
    for hdu in hdul:
        if hdu.data is not None:
            data = hdu.data
            break

hdul.close()

if data.ndim == 3:
    data = data[0]


vmin, vmax = np.percentile(data, [0, 100])

# cmap = plt.cm.gray.copy()
# cmap.set_bad(color="#9ECAE1")

masked = np.ma.masked_invalid(data)

vmin = np.nanmin(data)
vmax = np.nanmax(data)

plt.imshow(
    masked,
    cmap="gray",
    origin="lower",
    vmin=vmin,
    vmax=vmax
)

plt.axis("off")



plt.savefig(
    output,
    dpi=300,
    bbox_inches="tight",
    pad_inches=0
)

plt.close()

print("Done:", output)