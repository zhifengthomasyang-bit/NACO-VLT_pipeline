# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 14:12:48 2026

@author: yang.zhifeng
"""

import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve


# Richardson–Lucy deconvolution

def richardson_lucy(img, psf, iterations, eps=1e-8):
    img = np.nan_to_num(img).astype(np.float32)

    psf = np.nan_to_num(psf).astype(np.float32)
    psf = psf / (np.sum(psf) + eps)

    estimate = img.copy()

    psf_mirror = psf[::-1, ::-1]

    for i in range(iterations):

        conv = fftconvolve(estimate, psf, mode='same')


        relative_blur = img / (conv + eps)


        estimate *= fftconvolve(relative_blur, psf_mirror, mode='same')

    return estimate



img_path = r"D:\headers\alignment\calib 2004-01-12 Imaging epxt-0.4 OP3 FULL OP4 J HS DR S13\2004-01-13T03_01_51.396 - 2004-01-13T03_05_51.730.fits"
psf_path = r"D:\headers\PSF\PSF 2004-01-12 STD S165-E_9136_VIS expt-1.0 OP3 FULL OP4 J HS DR S13\PSF 2004-01-12 STD S165-E_9136_VIS expt-1.0 OP3 FULL OP4 J HS DR S13.fits"

img = fits.getdata(img_path)
psf = fits.getdata(psf_path)


img = np.nan_to_num(img)

img = img - np.median(img)
img[img < 0] = 0

from scipy.ndimage import shift

py, px = np.unravel_index(np.argmax(psf), psf.shape)
cy, cx = np.array(psf.shape) // 2

psf = shift(psf, (cy - py, cx - px))


result = richardson_lucy(img, psf, iterations=10)


fits.writeto(r"D:\headers\rl_result10 mean.fits", result, overwrite=True)

print("RL deconvolution done")