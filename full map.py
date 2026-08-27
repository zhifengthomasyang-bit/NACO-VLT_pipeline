# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 20:47:06 2026

@author: yang.zhifeng
"""

import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


folder = r"D:\headers\full map\2005\L_prime"

files = glob.glob(os.path.join(folder, "*.npz"))

print("number of maps:", len(files))

maps = []

for f in files:

    data = np.load(f)

    maps.append(
        {
            "file": f,
            "img": data["img"],
            "lon": data["lon"],
            "lat": data["lat"]
        }
    )


scale = np.ones(len(maps))


ref_id = 1
target_id = 2


ref = maps[ref_id]
target = maps[target_id]


img_ref = ref["img"]
img_target = target["img"]

overlap = (
    np.isfinite(img_ref)
    &
    np.isfinite(img_target)
)


print(
    "overlap pixels:",
    np.sum(overlap)
)



if np.sum(overlap) > 0:


    ratio = (
        img_ref[overlap]
        /
        img_target[overlap]
    )


    ratio = ratio[
        np.isfinite(ratio)
    ]


    ratio = ratio[
        (ratio > np.percentile(ratio,1))
        &
        (ratio < np.percentile(ratio,99))
    ]


    s = np.median(ratio)


    print(
        "scale factor:",
        s
    )


    scale[target_id] = s


else:

    print("No overlap found")

print("Final scales:")
for i,s in enumerate(scale):

    print(
        os.path.basename(files[i]),
        s
    )

d0 = np.load(files[0])

dlon = d0["lon"][0,1]-d0["lon"][0,0]
dlat = d0["lat"][1,0]-d0["lat"][0,0]

lon_start = 10

lon_global = np.arange(
    lon_start,
    lon_start+360+dlon,
    dlon
)

lat_global = np.arange(
    -90,
    90+dlat,
    dlat
)

sum_map = np.zeros(
    (len(lat_global),
     len(lon_global))
)

count_map = np.zeros_like(sum_map)


for i,m in enumerate(maps):


    print(
        "processing:",
        os.path.basename(m["file"])
    )


    img = m["img"].copy()
    lon = m["lon"].copy()
    lat = m["lat"]


    img *= scale[i]


    lon[lon < lon_start] += 360



    valid = np.isfinite(img)


    img = img[valid]
    lon = lon[valid]
    lat = lat[valid]



    lon_index = np.round(
        (lon-lon_start)/dlon
    ).astype(int)


    lat_index = np.round(
        (lat+90)/dlat
    ).astype(int)



    inside = (
        (lon_index>=0)
        &
        (lon_index<len(lon_global))
        &
        (lat_index>=0)
        &
        (lat_index<len(lat_global))
    )


    lon_index = lon_index[inside]
    lat_index = lat_index[inside]
    img = img[inside]



    np.add.at(
        sum_map,
        (lat_index,lon_index),
        img
    )


    np.add.at(
        count_map,
        (lat_index,lon_index),
        1
    )


merged = np.full_like(
    sum_map,
    np.nan
)


mask = count_map>0


merged[mask] = (
    sum_map[mask]
    /
    count_map[mask]
)



print(
    "coverage:",
    np.sum(mask)/mask.size
)


vmin=np.nanpercentile(
    merged,2
)

vmax=np.nanpercentile(
    merged,98
)


def format_longitude(x, pos):
    return f"{int(x % 360)}°"

plt.figure(figsize=(14,7))


plt.imshow(
    merged,
    origin="lower",
    extent=[
        lon_start,
        lon_start+360,
        -90,
        90
    ],
    cmap="gray",
    aspect="auto",
    vmin=vmin,
    vmax=vmax
)

ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_longitude))

plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")

plt.colorbar(
    label="Reflectance"
)

plt.title(
    "Ceres global map"
)

plt.show()