# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:47:11 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import os
from collections import Counter
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

path = r"D:\headers\catalogue\2004-01-15\FLAT,LAMP\classify expt\expt-5.0 MASK Slit_86mas OP3 FULL OP4 Grism4 OP6 SJ HS DR S54"
file_name = "FLAT expt-5.0 OP3 FULL OP4 Grism4 OP6 SJ HS DR S54.fits"
file_path = os.path.join(path, file_name)
cut1_deriv = 2.5
cut2_deriv = -2.5
cut_left = 10
cut_right = 5 
cut_bottum = 0
cut_up = 40

flat_data = []

with fits.open(file_path) as hdul:
    flat_data = hdul[0].data.copy()
    
flat_data = np.array(flat_data)

top_idx = np.argpartition(flat_data.ravel(), -100)[-100:]

rows, cols = np.unravel_index(top_idx, flat_data.shape)

median_row = np.median(rows)
median_col = np.median(cols)

median_col = int(median_col)
col_data = flat_data[:, median_col]
top_idx_col = np.argpartition(col_data, -50)[-50:]

std_row = np.median(top_idx_col)
std_row = int(std_row)

rows_data = flat_data[[std_row-20, std_row-10, std_row, std_row+10, std_row+20], :]
row_means = np.mean(rows_data, axis=1)


left_results = []
right_results = []
for i, row in enumerate(rows_data):

    cols = np.where(row < row_means[i])[0]

    left = []
    for c in cols:
        if c == len(left):
            left.append(c)
        else:
            break

    right = []
    for c in cols[::-1]:
        if c == 1023 - len(right):
            right.append(c)
        else:
            break

    if len(left) > 0:
        left_results.append(left[-1])

    if len(right) > 0:
        right_results.append(right[-1])


left_mode = Counter(left_results).most_common(1)[0][0]
right_mode = Counter(right_results).most_common(1)[0][0]

flat_masked = flat_data.copy()

flat_masked[:, 0:left_mode + cut_left] = np.nan
flat_masked[:, right_mode - cut_right:] = np.nan

col_cor = int((left_mode + right_mode) / 2 )

cols_cor = flat_data[: , [col_cor - 200 , col_cor - 100 , col_cor , col_cor + 100 , col_cor + 200 ]]

smooth = savgol_filter(cols_cor, window_length=21, polyorder=3, axis=0)

cols_deriv = np.gradient(smooth , axis=0)

row_index = np.arange(cols_cor.shape[0])

cut_positions = []

for i in range(cols_deriv.shape[1]):

    deriv_col = cols_deriv[:, i]

    idx = np.where(deriv_col > cut1_deriv)[0]

    if len(idx) > 0:
        cut_positions.append(int(idx[0]))
    else:
        cut_positions.append(None)

cut_row = int(np.median(cut_positions))

flat_masked[0:cut_row + cut_bottum, :] = np.nan


cut_positions2 = []

for i in range(cols_deriv.shape[1]):
    deriv_col = cols_deriv[:, i]
    
    idx = np.where(deriv_col < cut2_deriv)[0]    

    if len(idx) > 0:
        cut_positions2.append(int(idx[-1]))
    else:


        cut_positions2.append(cols_deriv.shape[0] - 1)



cut_row2 = int(np.median(cut_positions2))

flat_masked[cut_row2 - cut_up:, :] = np.nan

valid_rows_mask = ~np.isnan(flat_masked).all(axis=1)

row_medians = np.full(flat_masked.shape[0], np.nan) 
row_medians[valid_rows_mask] = np.nanmedian(flat_masked[valid_rows_mask], axis=1)

row_medians[row_medians == 0] = 1

flat_normalized = flat_masked / row_medians[:, np.newaxis]
flat_normalized = flat_normalized.astype(np.float32)

save_name = f"SPEC {file_name}"

save_path = os.path.join(path, save_name)

hdu = fits.PrimaryHDU(flat_normalized)

with fits.open(file_path) as hdul:
    hdu.header = hdul[0].header

hdu.writeto(save_path, overwrite=True)

