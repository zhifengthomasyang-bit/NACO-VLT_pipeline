# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:41:55 2026

@author: yang.zhifeng
"""

import os
import shutil
import csv
from astropy.io import fits

input_dir = r"D:\headers\catalogue\2005-05-09\STD\classify expt\expt-1.0"
output_dir = r"D:\headers\catalogue\2005-05-09\STD\classify expt\expt-1.0"

os.makedirs(output_dir, exist_ok=True)

results = []

def valid_filter(x):
    return x is not None and str(x).strip().lower() != "empty"

for filename in os.listdir(input_dir):
    if not filename.lower().endswith(".fits"):
        continue

    filepath = os.path.join(input_dir, filename)
    
    try:
        with fits.open(filepath) as hdul:
            header = hdul[0].header

            OP4 = header.get("HIERARCH ESO INS OPTI4 NAME")
            OP5 = header.get("HIERARCH ESO INS OPTI5 NAME")
            OP6 = header.get("HIERARCH ESO INS OPTI6 NAME")

            if valid_filter(OP4):
                op_name = "OP4"
                filt = OP4
            elif valid_filter(OP5):
                op_name = "OP5"
                filt = OP5
            elif valid_filter(OP6):
                op_name = "OP6"
                filt = OP6
            else:
                op_name = "OP_Unknow"
                filt = 'Unknow'
            
            filt = str(filt).strip()
                   
    except Exception:
        filt = "UNKNOWN"
    folder = f"OP3 FULL {op_name} {filt}"     
        
    dest_dir = os.path.join(output_dir, folder)
    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy2(filepath, dest_dir)