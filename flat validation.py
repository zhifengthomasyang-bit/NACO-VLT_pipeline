# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 11:58:43 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

with fits.open(r'D:\headers\master\FLAT,SKY\2004-01-12 MASK FLM_27\expt-3.5 OP3 FULL OP6 Ks HS DR S27.fits') as hdul:
    data = hdul[0].data
    

    p3 = np.nanpercentile(data, 3)   
    p99_9 = np.nanpercentile(data, 99.9) 
    
    print(f"90% data between: {p3:.4f} and {p99_9:.4f}")

valid_data = data[~np.isnan(data)]
p3 = np.nanpercentile(valid_data, 3)
p99_9 = np.nanpercentile(valid_data, 99.9)
subset = valid_data[(valid_data >= p3) & (valid_data <= p99_9)]


median_val = np.median(subset)
std_val = np.std(subset)


plt.figure(figsize=(10, 6))

count, bins, ignored = plt.hist(subset, bins=100, density=True, alpha=0.6, color='skyblue', label='LAMP Flat')


mode_val = bins[np.argmax(count)]
plt.axvline(median_val, color='red', linestyle='--', label=f'Median: {median_val:.2f}')
plt.axvline(mode_val, color='green', linestyle=':', label=f'Mode (Peak): {mode_val:.2f}')

plt.title(f'K-band LAMP Flat Probability Distribution (2.5%-99.9%)\nRange: [{p3:.2f}, {p99_9:.2f}]')
plt.xlabel('Data Units (DU)')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(axis='y', alpha=0.3)
plt.show()

peak_density_hist = np.max(count) 
print(f"max value: {peak_density_hist:.6f}")
print(f"Median: {median_val:.2f}")
print(f"Mode: {mode_val:.2f}")
print(f"Sigma: {std_val:.2f}")