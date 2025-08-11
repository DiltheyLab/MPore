#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:36:26 2024

@author: azlannisar
"""

import os
import pandas as pd

# List of barcodes
barcodes = ["barcode01","barcode02", "barcode03"]

# Base directories for the references and output
reference_base_dir = "/Users/azlannisar/Desktop/mount/MRSA_Data/Reference/"
output_base_dir = "/Users/azlannisar/Desktop/mount/MRSA_Data/Output_MRSA/"

# List to hold the data
data = []

# Loop through the barcodes and generate the paths
for barcode in barcodes:
    file_name = barcode  
    reference_path = os.path.join(reference_base_dir, barcode)
    pod5_path = os.path.join(output_base_dir, barcode)
    data.append([file_name,reference_path, pod5_path])

# Create a DataFrame
df = pd.DataFrame(data, columns=["File_name","Reference_path", "pod5_path"])

# Save the DataFrame to a CSV file
output_csv_path = "/Users/azlannisar/Snakemake_Output_final/MRSA_Data.csv"
df.to_csv(output_csv_path, index=False)

print(f"CSV file created at: {output_csv_path}")
