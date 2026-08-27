# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:51:55 2026

@author: yang.zhifeng
"""

import os
import re
import pandas as pd

folder_path = r"D:\headers\catalogue\2005-05-09\DARK expt-0.2 HWD Un"
all_data = []

for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)

        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        data = {}

        matches = re.findall(r'^([A-Z0-9 _]+):\s*(.*)$', content, re.MULTILINE)
        for key, value in matches:
            key = key.strip().replace(" ", "_")
            value = value.strip()

            if re.fullmatch(r'OPTI[1-7]', key, re.IGNORECASE):
                continue

            data[key] = value

        match = re.search(r'^\s*OPTI1:?\s*(.*)$', content, re.MULTILINE | re.IGNORECASE)
        if match:
            data['MASK'] = match.group(1).strip()

        for i in range(2, 7):
            match = re.search(rf'^\s*OPTI{i}:?\s*(.*)$', content, re.MULTILINE | re.IGNORECASE)
            if match:
                data[f'OPTI{i}'] = match.group(1).strip()

        match = re.search(r'^\s*OPTI7:?\s*(.*)$', content, re.MULTILINE | re.IGNORECASE)
        if match:
            data['CAMERA'] = match.group(1).strip()

        data["FILENAME"] = filename

        all_data.append(data)

priority_cols = ["FILENAME", "MASK", "OPTI2", "OPTI3", "OPTI4", "OPTI5", "OPTI6", "CAMERA"]


all_keys = []
for row in all_data:
    for key in row.keys():
        if key not in all_keys:
            all_keys.append(key)

all_keys = [col for col in priority_cols if col in all_keys] + [col for col in all_keys if col not in priority_cols]

df = pd.DataFrame(all_data, columns=all_keys)

df = df.fillna("").astype(str)

max_width = {}
for col in df.columns:
    max_len = max(df[col].map(len).max(), len(col))
    max_width[col] = max_len + 4   

for col in df.columns:
    width = max_width[col]
    df[col] = df[col].apply(lambda x: x.center(width))


new_columns = {col: col.center(max_width[col]) for col in df.columns}
df.rename(columns=new_columns, inplace=True)

output_path = "output.csv"
df.to_csv(output_path, index=False, encoding="utf-8-sig")
