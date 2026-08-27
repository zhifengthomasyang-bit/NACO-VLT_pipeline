# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:24:16 2026

@author: yang.zhifeng
"""

from PIL import Image

img1 = Image.open(r"D:\headers\Test\orthographic projection.png").convert("RGB")
img2 = Image.open(r"D:\headers\Test\orthographic projection 1.png").convert("RGB")

img2 = img2.resize(img1.size)

img1.save(
    "compare.gif",
    save_all=True,
    append_images=[img2],
    duration=1000,
    loop=0,
    optimize=False
)

