# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:09:56 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
import glob

dark_files = (#glob.glob(r"D:\headers\catalogue\2004-01-12\FLAT,LAMP\classify expt\expt-3.0\OP4 J HS DR S27\DARK\*.fits")
              #+glob.glob(r"D:\headers\catalogue\2004-01-12\FLAT,LAMP\classify expt\expt-3.0\OP6 H HS DR S27\DARK\*.fits")
              glob.glob(r"D:\headers\catalogue\2005-05-09\DARK expt-0.2 HWD Un\*.fits")
              )
#dark_files = [
    #r"D:\headers\catalogue\2004-01-12\FLAT,LAMP\classify expt\expt-3.0\OP5 NB_2.17 HS FN S27\DARK\NACO.2004-01-13T10%3A36%3A03.181.fits",
    #r"D:\headers\catalogue\2004-01-12\FLAT,LAMP\classify expt\expt-3.0\OP5 NB_2.17 HS FN S27\DARK\NACO.2004-01-13T10%3A38%3A13.086.fits",
    #r"D:\headers\catalogue\2004-01-12\FLAT,LAMP\classify expt\expt-3.0\OP5 NB_2.17 HS FN S27\DARK\NACO.2004-01-13T10%3A40%3A24.072.fits"]

dark_stack = []

for file in dark_files:
    with fits.open(file) as hdul:
        data = hdul[0].data.astype(np.float32)
        dark_stack.append(data)

dark_stack = np.array(dark_stack)


master_dark = np.median(dark_stack, axis=0)


hdu = fits.PrimaryHDU(master_dark)
hdu.writeto("expt-0.2 HWD Un L27.fits", overwrite=True)

