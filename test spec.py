# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:14:20 2026

@author: yang.zhifeng
"""

import matplotlib.pyplot as plt
from astropy.io import fits
import os


file_path = r"D:\headers\catalogue\2004-01-12\WAVE,LAMP\SJ\FLAT SJ.fits"

with fits.open(file_path) as hdul:
    data = hdul[0].data
    row_512 = data[628, :]
    col_512 = data[:, 511]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].plot(row_512, color='blue', linewidth=1)
    axes[0].set_title("ADU Distribution: Row 512")
    axes[0].set_xlabel("Column Index")
    axes[0].set_ylabel("ADU Value")
    
    
    axes[1].plot(col_512, color='red', linewidth=1)
    axes[1].set_title("ADU Distribution: Column 512")
    axes[1].set_xlabel("Row Index")
    axes[1].set_ylabel("ADU Value")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()