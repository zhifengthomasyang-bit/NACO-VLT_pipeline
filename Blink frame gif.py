# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:34:27 2026

@author: yang.zhifeng
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

files = [r"D:\headers\catalogue\2004-01-12\OBJECT\classify expt\expt-1.0 MASK Slit_86mas OP3 FULL OP4 Grism2 OP6 SK HS DR S54\A\NACO.2004-01-13T04%3A49%3A17.947.fits",
    r"D:\headers\catalogue\2004-01-12\OBJECT\classify expt\expt-1.0 MASK Slit_86mas OP3 FULL OP4 Grism2 OP6 SK HS DR S54\B\NACO.2004-01-13T04%3A50%3A46.123.fits",
    r"D:\headers\calib\2004-01-12\Spec\expt-1.0 MASK Slit_86mas OP3 FULL OP4 Grism2 OP6 SK HS DR S54\Calib NACO.2004-01-13T04_49_17.947.fits - NACO.2004-01-13T04_50_46.123.fits",
]


images = [
    np.squeeze(fits.getdata(f)).astype(float)
    for f in files
]

allpix = np.concatenate([im.ravel() for im in images])
vmin, vmax = np.nanpercentile(allpix, [0, 100])

fig, ax = plt.subplots(figsize=(8, 8))

im = ax.imshow(
    images[0],
    origin="lower",
    cmap="gray",
    vmin=vmin,
    vmax=vmax
)

ax.set_axis_off()

def update(i):
    im.set_data(images[i])
    return [im]

ani = FuncAnimation(
    fig,
    update,
    frames=len(images),
    interval=250,
    blit=True,
    repeat=True
)

ani.save("A B A-B.gif", writer="pillow", fps=1)