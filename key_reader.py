# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:42:51 2026

@author: yang.zhifeng
"""
import contextlib
from astropy.io import fits
import matplotlib.pyplot as plt
import os

def list_all_filters(header):
    print("=== All FILTER elements ===")
    
    for i in range(1, 10):
        type_key = f'HIERARCH ESO INS OPTI{i} TYPE'
        name_key = f'HIERARCH ESO INS OPTI{i} NAME'
        
        if type_key in header and str(header[type_key]).strip() == 'FILTER':
            name = str(header.get(name_key, 'Unknown')).strip()
            print(f"OPTI{i}: {name}")

def find_camera(header):
    print("=== Camera ===")
    
    for i in range(1, 10):
        type_key = f'HIERARCH ESO INS OPTI{i} TYPE'
        name_key = f'HIERARCH ESO INS OPTI{i} NAME'
        
        if type_key in header and str(header[type_key]).strip() == 'OBJECTIVE':
            name = str(header.get(name_key, 'Unknown')).strip()
            print(f"OPTI{i}: {name}")
        
def find_mask(header):
    print("=== Mask ===")
    
    for i in range(1, 10):
        type_key = f'HIERARCH ESO INS OPTI{i} TYPE'
        name_key = f'HIERARCH ESO INS OPTI{i} NAME'
        
        if type_key in header and str(header[type_key]).strip() == 'MASK':
            name = str(header.get(name_key, 'Unknown')).strip()
            print(f"OPTI{i}: {name}")
    
        
def read_fits(filename):
    hdul = fits.open(filename)
    header = hdul[0].header
    keys_of_interest = ['EXPTIME', 'OBJECT', 'DATE-OBS', 'HIERARCH ESO DPR TYPE', 'HIERARCH ESO DPR CATG','HIERARCH ESO DET MODE NAME','HIERARCH ESO DET NCORRS NAME']
    display_name = os.path.basename(filename).replace('%3A', ':')
    print(f"=== key words ({display_name}) ===")
    for key in keys_of_interest:
        if key in header:
            print(f"{key}: {header[key]}")
        else:
            print(f"{key}: unexist")
    list_all_filters(header)
    find_camera(header)
    find_mask(header)

    hdul.close()
    
    return header

filename = r'D:\headers\catalogue\2005-05-09\NACO.2005-05-10T13_03_29.057.fits'
output = r'D:\headers'

read_fits(filename)

if not os.path.exists(output):
    os.makedirs(output)
txt_file = os.path.join(output, os.path.splitext(os.path.basename(filename).replace('%3A', '_'))[0] + '.txt')


with open(txt_file, 'w') as f:
    with contextlib.redirect_stdout(f):
        read_fits(filename)