# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 13:36:15 2026

@author: yang.zhifeng
"""


import numpy as np
from astropy.io import fits
from scipy.signal import wiener, fftconvolve

img = fits.getdata(r"D:\headers\alignment\calib 2004-01-12 Imaging epxt-0.4 OP3 FULL OP4 J HS DR S13\2004-01-13T03_01_51.396 - 2004-01-13T03_05_51.730.fits").astype(np.float32)
psf = fits.getdata(r"D:\headers\PSF\PSF 2004-01-12 STD S165-E_9136_VIS expt-0.5 OP3 FULL OP4 J HS DR S27\PSF 2004-01-12 STD S165-E_9136_VIS expt-0.5 OP3 FULL OP4 J HS DR S27.fits").astype(np.float32)

img = np.nan_to_num(img)
psf = psf / np.sum(psf)

# Wiener 反卷积（简化版）
def wiener_deconv(img, psf, K=0.01):
    psf_ft = np.fft.fft2(psf, s=img.shape)
    img_ft = np.fft.fft2(img)

    psf_conj = np.conj(psf_ft)
    result_ft = img_ft * psf_conj / (np.abs(psf_ft)**2 + K)

    return np.abs(np.fft.ifft2(result_ft))

result = wiener_deconv(img, psf)

fits.writeto("D:\\headers\\wiener_result.fits", result, overwrite=True)