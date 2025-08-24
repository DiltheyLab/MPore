#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 17:20:47 2024

@author: azlannisar
"""

import pandas as pd
import os

# Get input CSV file and output directory from environment variables
input_csv = os.getenv("INPUT_CSV")
output_dir = os.getenv("OUTPUT_DIR")

# Read the input CSV file
data = pd.read_csv(input_csv)

# Add a column called 'Bed_path' and populate it with values
data['Bam_data'] = output_dir + '/' + data['File_name'] + '_sorted.bam'
data['Bed_data'] = output_dir + '/' + data['File_name'] + '.bed'

input_base = os.path.splitext(os.path.basename(input_csv))[0]
updated_csv = os.path.join(output_dir, f"{input_base}_updated.csv")


data.to_csv(updated_csv, index=False)