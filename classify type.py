# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 09:29:43 2026

@author: yang.zhifeng
"""

import pandas as pd

df = pd.read_csv(
    "D:\\headers\\catalogue\\2005-05-09\\catalogue-2005-05-09.csv",
    dtype=str)

print(df["  HIERARCH_ESO_DPR_TYPE  "].unique())

for t, group in df.groupby("  HIERARCH_ESO_DPR_TYPE  "):
    group.to_csv(f"{t}.csv", index=False)
    
    
    