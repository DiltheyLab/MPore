#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 17:56:05 2024

@author: azlannisar
"""

import os
import glob
import matplotlib.pyplot as plt
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
    
    # Get unique motifs from the first DataFrame
    unique_motifs = dataframes_list[0]['Motif'].unique()
    
    # Check if there are files in matched_files_5mC
    if matched_files_5mC:
        # Create a dictionary to hold scores for each motif
        motif_scores = {}
    
        # Iterate over each motif and collect corresponding scores
        for motif in unique_motifs:
            scores = []
            for df in dataframes_list:
                df_motif = df[df['Motif'] == motif]
                scores.append(df_motif['Score1'].values)
            motif_scores[motif] = scores
    
        # Plotting one boxplot for each unique motif across all DataFrames
        for motif, scores in motif_scores.items():
            if 'C' in motif:
                clean_motif = motif.replace('"', '')
                matching_rows = enzyme_presence_df[enzyme_presence_df['Motif'].str.contains(clean_motif.split(',')[0], na=False)]
                fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the figure size as needed
                box = ax.boxplot(scores, labels=x_values, showfliers=False)
                ax.set_xlabel('')
                ax.set_ylabel('Score')
                ax.set_ylim(0, 110)  # Set the Y-axis limits from 0 to 100
                ax.set_title(f'5mC Boxplot for {motif}')
                plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability 
                
                ax.axhline(y=25, color='magenta', linestyle='--', linewidth=1)
    
                if not matching_rows.empty:
                    # Create a DataFrame for the enzyme presence with e-values
                    table_data = matching_rows.values
                    col_labels = matching_rows.columns
                    #cell_text = [[f"{value}" for value in row] for row in table_data]
                    cell_text = []
                    for row in table_data:
                        formatted_row = []
                        for value in row:
                            if isinstance(value, (float, np.float64)):
                                # Format float values to retain significant digits
                                formatted_value = '{:.10g}'.format(value)
                            else:
                                formatted_value = str(value)
                            formatted_row.append(formatted_value)
                        cell_text.append(formatted_row)
    
                    # Add the table to the figure
                    table = plt.table(cellText=cell_text, colLabels=col_labels, cellLoc='center', loc='bottom', bbox=[-0.1, -1.7, 1.1, 1.1])
                    table.auto_set_font_size(False)
                    table.set_fontsize(8)
                    table.scale(1, 1.5)
    
                    # Adjust layout to make space for the table
                    plt.subplots_adjust(left=0.2, bottom=0.6)  # Adjust bottom to make space for the table
                else:
                # Adjust layout when no table is added
                    plt.subplots_adjust(left=0.2, bottom=0.2)  # Adjust bottom to make space for the table
                #plt.show()
                plt.savefig(os.path.join(args.output_dir, f'{motif}_5mC_boxplot.png'))
                plt.close()
        
    else:
        print("No 5mC Files found.")       

if __name__ == "__main__":
    main()