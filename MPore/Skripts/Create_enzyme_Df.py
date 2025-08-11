#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 17:32:01 2024

@author: azlannisar
"""

import os
import glob
import pandas as pd
import argparse
import numpy as np

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Mtase_presence_Plot1')
    parser.add_argument('MTASE_FILE', type=str, help='Path_to_Mtase_file_path')
    parser.add_argument('Input_dir', type=str, help='Path to list dir')
    parser.add_argument('output_dir', type=str, help='Output directory path')
    return parser.parse_args()

def main():
    args = parse_arguments()
    directory = args.Input_dir
    #directory= "/Users/azlannisar/Snakemake_Test"
    # Load 5mC Files
    matched_files_5mC = glob.glob(os.path.join(directory, '*_detailed_*_5mC.csv'))
    dataframes_list = []
    for file_path in matched_files_5mC:
        # Read CSV file into a DataFrame
        df = pd.read_csv(file_path)
        df = df[df['Score1'] != '/']
        df['Score1'] = pd.to_numeric(df['Score1'], errors='coerce')
        # Append DataFrame to the list
        dataframes_list.append(df)
    
    # Read Mtase_presence DataFrame
    mtase_presence_path = args.MTASE_FILE
    #mtase_presence_path = "/Users/azlannisar/Mtase_presence_e_25_values.csv"
    Mtase_presence = pd.read_csv(mtase_presence_path, sep=',', index_col=False)
    
    tsv_file_path = "./Type2_Datframe.tsv"
    tsv_df = pd.read_csv(tsv_file_path, sep='\t')
    x_values = [os.path.basename(file_path).split('_detailed_')[1].split('_')[0] for file_path in matched_files_5mC]
    
    enzyme_presence_df = pd.DataFrame(columns=['Enzyme', 'Motif'] + x_values)
    
    
    
    # Populate the DataFrame with e-values and add Mod_Position to Motif column
    for enzyme in Mtase_presence.columns[1:]:  # Skip the 'Isolate' column
        row = [enzyme, np.nan]  # Initialize Motif column with NaN
        for x_label in x_values:
            isolate_df = Mtase_presence[Mtase_presence['Isolates'] == x_label]
            if not isolate_df.empty:
                row.append(isolate_df[enzyme].values[0])
            else:
                row.append('NA')
        enzyme_presence_df.loc[len(enzyme_presence_df)] = row
    
    enzyme_presence_df = enzyme_presence_df.merge(tsv_df[['Enzyme', 'Mod_Position']], on='Enzyme', how='left')
    enzyme_presence_df['Motif'] = enzyme_presence_df['Mod_Position']
    enzyme_presence_df.drop(columns=['Mod_Position'], inplace=True)  # Remove Mod_Position column after merging
    
    enzyme_presence_output_path = os.path.join(args.output_dir, "enzyme_presence_df.csv")
    enzyme_presence_df.to_csv(enzyme_presence_output_path, index=False)

if __name__ == "__main__":
    main()