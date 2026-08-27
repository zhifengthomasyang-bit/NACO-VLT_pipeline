# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 10:54:17 2026

@author: yang.zhifeng
"""

import contextlib
from astropy.io import fits
import os
import key_reader  

input_folder = 'D:\\ceres\\2004-01-15b\\'
output_folder = 'D:\\headers\\2004-01-15b\\'



if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 遍历整个文件夹里的所有 fits 文件
for file in os.listdir(input_folder):
    if file.lower().endswith('.fits'):
        filename = os.path.join(input_folder, file)
        txt_file = os.path.join(output_folder, os.path.splitext(file.replace('%3A', '_'))[0] + '.txt')
        
        # 打开 txt 文件，把 print 内容写入
        with open(txt_file, 'w') as f:
            with contextlib.redirect_stdout(f):
                key_reader.read_fits(filename)
        
        print(f"Processed: {file}")