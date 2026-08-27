
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:09:33 2026

@author: yang.zhifeng
"""

import pandas as pd
import os
from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clip
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from numpy.lib.stride_tricks import sliding_window_view
from matplotlib.patches import Circle
from scipy.interpolate import griddata

def csv_to_dict(path):
    df = pd.read_csv(path, nrows=2, header=None, sep=';')

    indices = df.iloc[0]
    values = df.iloc[1]

    mapping = dict(zip(indices, values))

    out = {}

    for k, v in mapping.items():
        key = (
            str(k)
            .replace(' ', '_')
            .replace('(', '')
            .replace(')', '')
            .replace('/', '_')
            .replace('.', '')
            .replace('-', '_')
        )
        
        try:
            v = float(v)
        except:
            pass
        out[key] = v
    return out

variables = {}

path = r"D:\headers\alignment\calib 2005-05-08 Imaging Ceres8_K expt-2.0 OP3 FULL OP6 Ks HD DR S13"
variables.update(csv_to_dict(os.path.join(path, "Ephemcc 2005-05-09T02_27_03.635.csv")))
variables.update(csv_to_dict(os.path.join(path, "Ephemph 2005-05-09T02_27_03.635.csv")))
fits_path = os.path.join(path, "2005-05-09T02_20_20.703 - 2005-05-09T02_34_27.306.fits")

path2 = r"D:\headers\Test"

with fits.open(fits_path) as hdul:
    img = np.squeeze(hdul[0].data).astype(float)
    

pixel_scale = 0.013260

SEP_long_deg = float(variables['SEP_long_deg'])
SEP_lat_deg = float(variables['SEP_lat_deg'])
SSP_long_deg = float(variables['SSP_long_deg'])
SSP_lat_deg = float(variables['SSP_lat_deg'])
NP_deg = float(variables['NP_deg'])
dPole_arcsec = float(variables['dPole_arcsec'])
V_Mag_mag = float(variables['V_Mag_mag'])
Phase_deg = float(variables['Phase_deg'])
App_Radius_arcsec = float(variables['App_Radius_arcsec'])
Dobs_au = float(variables['Dobs_au'])
Dh_au = float(variables['Dh_au'])
PAQ_deg = float(variables['PAQ_deg'])
Q_mas = float(variables['Q_mas'])
RA_h = float(variables['RA_h'])
DEC_deg = float(variables['DEC_deg'])
VMag_mag = float(variables['VMag_mag'])
Elong_deg = float(variables['Elong_deg'])
dRAcosDEC_arcsec_min = float(variables['dRAcosDEC_arcsec_min'])
dDEC_arcsec_min = float(variables['dDEC_arcsec_min'])
RV_km_s = float(variables['RV_km_s'])

AU = 1.495978707e8

R_pixel = App_Radius_arcsec / pixel_scale  

y_size, x_size = img.shape
x0, y0 = x_size / 2.0, y_size / 2.0

R_pixel_prime = (np.arctan(483.5 / (Dobs_au * AU)) * 648000 / np.pi) / pixel_scale
coeff = 446.0 / 483.5

lambda_0 = np.radians(SEP_long_deg)
phi_0 = np.radians(SEP_lat_deg)
P =  np.radians(NP_deg)

def geo_to_pixel(lon_deg, lat_deg):

    lon_rad = np.radians(lon_deg)
    lat_rad = np.radians(lat_deg)
    

    X_body = np.cos(lat_rad) * np.sin(lon_rad - lambda_0)
    Z_body = np.cos(lat_rad) * np.cos(lon_rad - lambda_0)
    Y_body = np.sin(lat_rad) * coeff
    
    x_prime = X_body
    y_prime = Y_body * np.cos(phi_0) - Z_body * np.sin(phi_0)
    z_prime = Y_body * np.sin(phi_0) + Z_body * np.cos(phi_0)
    

    if np.any(z_prime <= 0):
        pass 
        

    x_disk = x_prime * np.cos(P) - y_prime * np.sin(P)
    y_disk = x_prime * np.sin(P) + y_prime * np.cos(P)
    

    px = x0 + x_disk * R_pixel_prime
    py = y0 + y_disk * R_pixel_prime
        
    if isinstance(px, np.ndarray):
        visible = z_prime > 0
        px = px.astype(float) 
        py = py.astype(float)
        px[~visible] = np.nan
        py[~visible] = np.nan
    else:
        if z_prime <= 0:
            return np.nan, np.nan
            
    return px, py


ssp_x, ssp_y = geo_to_pixel(SSP_long_deg, SSP_lat_deg)


plt.figure(figsize=(9, 9), dpi=100)


v_min = np.min(img)
v_max = np.max(img)
plt.imshow(img, cmap='gray', origin='lower', vmin=v_min, vmax=v_max)


lat_lines = [-60, -30, 0, 30, 60]
lon_samples = np.linspace(0, 360, 500) 

for lat in lat_lines:
    px, py = geo_to_pixel(lon_samples, np.full_like(lon_samples, lat))

    if lat == 0:
        plt.plot(px, py, 'r-', linewidth=1.5, alpha=0.6, label='Equator (0°)')
    else:
        plt.plot(px, py, 'r:', linewidth=1.0, alpha=0.4)


lon_lines = np.arange(0, 360, 30)
lat_samples = np.linspace(-90, 90, 500)

for lon in lon_lines:
    px, py = geo_to_pixel(np.full_like(lat_samples, lon), lat_samples)
    plt.plot(px, py, 'r:', linewidth=1.0, alpha=0.4)

lon_grid = np.linspace(-180, 180, 720)
lat_grid = np.linspace(-90, 90, 360)

LON, LAT = np.meshgrid(lon_grid, lat_grid)


# project surface grid
px, py = geo_to_pixel(
    LON,
    LAT
)


mask = np.isfinite(px) & np.isfinite(py)


x_vis = px[mask]
y_vis = py[mask]


from scipy.spatial import ConvexHull

points = np.column_stack((x_vis, y_vis))

hull = ConvexHull(points)

limb = points[hull.vertices]


plt.plot(
    np.append(limb[:,0], limb[0,0]),
    np.append(limb[:,1], limb[0,1]),
    'g--',
    linewidth=1.5,
    label='Ceres Ellipsoid Limb'
)
# theta = np.linspace(0, 2*np.pi, 200)
# plt.plot(x0 + R_pixel_prime*np.cos(theta), y0 + R_pixel_prime*np.sin(theta), 'g--', alpha=0.7, label='Ceres Limb')


plt.plot(x0, y0, 'bx', markersize=8, markeredgewidth=2, label=f'SEP (Sub-Earth): {SEP_long_deg}°, {SEP_lat_deg}°')


plt.plot(ssp_x, ssp_y, 'r*', markersize=14, markeredgecolor='white', label=f'SSP (Sub-Solar): {SSP_long_deg}°, {SSP_lat_deg}°')


plt.title("Ceres NACO Image with Orthographic Lat/Lon Grid", fontsize=12)
plt.legend(loc='upper right', fontsize=9)
plt.colorbar(label='Intensity (ADU)', fraction=0.046, pad=0.04)
plt.grid(False) 

plt.savefig(os.path.join(path2, "orthographic projection.png"))
plt.show()


yy, xx = np.meshgrid(np.arange(y_size), np.arange(x_size), indexing='ij')


disk_mask = (xx - x0)**2 + (yy - y0)**2 <= R_pixel**2

lon_test = np.linspace(SEP_long_deg - 120,
                       SEP_long_deg + 120,
                       180)
lon_test = lon_test % 360

lat_test = np.linspace(-90, 90, 90)

lon_grid, lat_grid = np.meshgrid(lon_test, lat_test)

px_model, py_model = geo_to_pixel(lon_grid, lat_grid)


px_f = px_model.flatten()
py_f = py_model.flatten()
lon_f = lon_grid.flatten()
lat_f = lat_grid.flatten()


lon_map = np.full(img.shape, np.nan)
lat_map = np.full(img.shape, np.nan)

for i in range(y_size):
    for j in range(x_size):

        if not disk_mask[i, j]:
            continue

        x = xx[i, j]
        y = yy[i, j]

        dist = (px_f - x)**2 + (py_f - y)**2
        idx = np.nanargmin(dist)

        lon_map[i, j] = lon_f[idx]
        lat_map[i, j] = lat_f[idx]


lon_valid = lon_map[disk_mask]
lat_valid = lat_map[disk_mask]
img_valid = img[disk_mask]

valid = ~np.isnan(lon_valid) & ~np.isnan(lat_valid)

lon_valid = lon_valid[valid]
lat_valid = lat_valid[valid]
img_valid = img_valid[valid]


lon_min, lon_max = np.nanpercentile(lon_valid, [2, 98])
lat_min, lat_max = np.nanpercentile(lat_valid, [2, 98])

print("Longitude range:", lon_min, lon_max)
print("Latitude range:", lat_min, lat_max)


lon_out = np.linspace(lon_min, lon_max, 500)
lat_out = np.linspace(lat_min, lat_max, 250)

lon2d, lat2d = np.meshgrid(lon_out, lat_out)


img_lonlat = griddata(
    (lon_valid, lat_valid),
    img_valid,
    (lon2d, lat2d),
    method='linear'
)


plt.figure(figsize=(12, 5))

plt.imshow(
    img_lonlat,
    extent=[lon_min, lon_max, lat_min, lat_max],
    origin='lower',
    cmap='gray',
    aspect='auto'
)


plt.plot(SEP_long_deg, SEP_lat_deg,
            'bx', markersize=8, 
            label=f'SEP: {SEP_long_deg}°, {SEP_lat_deg}°')

plt.plot(SSP_long_deg,SSP_lat_deg,
            'r*', markersize=14,
            label=f'SSP: {SSP_long_deg}°, {SSP_lat_deg}°')
plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")
plt.title("Ceres Local Lat/Lon Unwrapped Map (Limb-constrained)")
plt.colorbar(label="Intensity (ADU)")

plt.legend(loc='upper right', fontsize=9)

plt.savefig(os.path.join(path2, "Unwrapped_Map_2-98.png"))
plt.show()


angle = 60.0
angle_deg = (90.0 - angle) * np.pi / 180.0
AU = 1.495978707e8
App_Radius_rad = App_Radius_arcsec * np.pi / 648000
Radius = np.tan(App_Radius_rad) * Dobs_au * AU
D_prime = - np.sin(angle_deg) * Radius  + np.sqrt(4 * (Dobs_au * AU) ** 2 - 4 * Radius ** 2 + (2 * Radius * np.sin(angle_deg))) / 2
cos_angle_eff = (Radius ** 2 + (Dobs_au * AU) ** 2 - D_prime ** 2) / (2 * Radius * (Dobs_au * AU))
angle_eff = np.arccos(cos_angle_eff)
angle_eff_deg = angle_eff * 180.0 / np.pi

deg_rad = np.pi / 180.0

phi = lat2d * deg_rad
lam = lon2d * deg_rad

phi0 = SEP_lat_deg * deg_rad
lam0 = SEP_long_deg * deg_rad

cos_theta = (
    np.sin(phi) * np.sin(phi0) +
    np.cos(phi) * np.cos(phi0) * np.cos(lam - lam0)
)

theta = np.arccos(cos_theta)
theta_deg = theta * 180.0 / np.pi

mask_60 = theta_deg < angle_eff_deg
img_masked = np.where(mask_60, img_lonlat, np.nan)

plt.figure(figsize=(12, 5))

plt.imshow(
    img_masked,
    extent=[lon_min, lon_max, lat_min, lat_max],
    origin='lower',
    cmap='gray',
    aspect='auto'
)

plt.plot(SEP_long_deg, SEP_lat_deg, 'bx', markersize=8, label='SEP')
plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")
plt.title(f"Ceres Map (emission angle < {angle}°)")
plt.colorbar(label="Intensity (ADU)")
plt.legend()

plt.savefig(os.path.join(path2, f"emission {angle}.png"))
plt.show()

phi1 = SSP_lat_deg * deg_rad
lam1 = SSP_long_deg * deg_rad

cos_theta1 = (
    np.sin(phi) * np.sin(phi1) +
    np.cos(phi) * np.cos(phi1) * np.cos(lam - lam1)
)
theta1 = np.arccos(cos_theta1)
theta_deg1 = theta1 * 180.0 / np.pi

mask_60 = theta_deg1 < angle_eff_deg
img_masked = np.where(mask_60, img_lonlat, np.nan)

plt.figure(figsize=(12, 5))

plt.imshow(
    img_masked,
    extent=[lon_min, lon_max, lat_min, lat_max],
    origin='lower',
    cmap='gray',
    aspect='auto'
)

plt.plot(SSP_long_deg, SSP_lat_deg, 'r*', markersize=14, label='SEP')
plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")
plt.title(f"Ceres Map (incidence angle < {angle}°)")
plt.colorbar(label="Intensity (ADU)")
plt.legend()

plt.savefig(os.path.join(path2, f"incidence {angle}.png"))
plt.show()