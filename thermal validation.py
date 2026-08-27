# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:36:47 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np

J_file = [r"D:\headers\master\DARK,LAMP\2004-01-12\expt-15.0 OP4 J HS DR S13.fits"]
J_dark = []

H_file = [r"D:\headers\master\DARK,LAMP\2004-01-12\expt-15.0 OP6 H HS DR S13.fits"]
H_dark =[]

Ks_file = [r"D:\headers\master\DARK,LAMP\2004-01-12\expt-15.0 OP6 Ks HS DR S13.fits"]
Ks_dark = []

for f in J_file:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    J_dark.append(data)
    
J_dark = np.array(J_dark)

for f in H_file:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    H_dark.append(data)
    
H_dark = np.array(H_dark)

for f in Ks_file:
    hdul = fits.open(f)
    data = hdul[0].data.astype(np.float32)
    Ks_dark.append(data)
    
Ks_dark = np.array(Ks_dark)

J_minus_H = J_dark - H_dark
J_minus_Ks = J_dark - Ks_dark

hdu = fits.PrimaryHDU(J_minus_H)
hdu.writeto("J_minus_H.fits", overwrite=True)

hdu = fits.PrimaryHDU(J_minus_Ks)
hdu.writeto("J_minus_Ks.fits", overwrite=True)