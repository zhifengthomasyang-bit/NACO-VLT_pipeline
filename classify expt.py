# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:49:56 2026

@author: yang.zhifeng
"""

import os
import shutil
import csv
from astropy.io import fits

input_dir = r"D:\headers\catalogue\2005-05-09\STD"
output_dir = r"D:\headers\catalogue\2005-05-09\STD\classify expt"

os.makedirs(output_dir, exist_ok=True)

results = []

for filename in os.listdir(input_dir):
    if not filename.lower().endswith(".fits"):
        continue

    filepath = os.path.join(input_dir, filename)

    try:
        with fits.open(filepath) as hdul:
            header = hdul[0].header

            exptime = header.get("EXPTIME")
            if exptime is None and len(hdul) > 1:
                exptime = hdul[1].header.get("EXPTIME")

    except Exception as e:
        print(f"failed: {filename}, error: {e}")
        exptime = None
        
    exptime_val = float(exptime)

    if exptime_val.is_integer():
        exptime_str = f"{exptime_val:.1f}"
    else:
        exptime_str = str(exptime_val)
    
    category = f"expt-{exptime_str}"

    category_dir = os.path.join(output_dir, category)
    os.makedirs(category_dir, exist_ok=True)

    dst_path = os.path.join(category_dir, filename)
    shutil.copy2(filepath, dst_path)

    results.append([filename, exptime, category, dst_path])

