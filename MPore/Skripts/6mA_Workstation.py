#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 21:19:17 2024

@author: azlannisar
"""

def extract_info(position):
    # Check if position is NaN
    if pd.isna(position):
        return [np.nan], [np.nan], [np.nan]
    else:
        # Convert position to string to ensure it is treated as a string
        position = str(position)
        
        positions = []
        methyl_types = []
        strands = []
        
        # Split positions by comma
        entries = position.split(',')
        
        for entry in entries:
            entry = entry.strip()  # Remove any leading/trailing whitespace
            
            if entry.startswith('-'):
                strand = '-'
                pos, methyl_type = entry[1:].split('(')
            else:
                strand = '+'
                pos, methyl_type = entry.split('(')
            
            pos = pos.strip()
            methyl_type = methyl_type.strip(')')
            methyl_type = methyl_type + 'mC' if methyl_type in ['5', '4'] else methyl_type + 'mA'
            
            positions.append(pos)
            methyl_types.append(methyl_type)
            strands.append(strand)
        
        return positions, methyl_types, strands
def generate_combinations(motif, methyl_type, enzyme, strand, known_position=None):
    results = []
    if methyl_type == '5mC':
        base = 'C'
    elif methyl_type == '6mA':
        base = 'A'
    elif methyl_type == '4mC':
        base = 'C'
    else:
        return results

    if known_position is None:
        for i, nucleotide in enumerate(motif):
            if nucleotide == base:
                unique_name = f"{enzyme}_{base}{i+1}"
                results.append([enzyme, unique_name, motif, methyl_type, strand, str(i+1)])
    else:
        unique_name = f"{enzyme}_{base}{known_position}"
        results.append([enzyme, unique_name, motif, methyl_type, strand, str(known_position)])
    
    return results
def process_data(df):
    result_rows = []

    for _, row in df.iterrows():
        enzyme = row['Enzyme']
        motif = row['Motif']
        positions = str(row['NewPosition']).split(',') if not pd.isna(row['NewPosition']) else [np.nan]
        methyl_types = str(row['MethylType']).split(',') if not pd.isna(row['MethylType']) else [np.nan]
        strands = str(row['Strand']).split(',') if not pd.isna(row['Strand']) else [np.nan]

        # Case 1: Klare Position
        if len(positions) == 1 and positions[0] != '?':
            result_rows.append([enzyme, enzyme, motif, methyl_types[0], strands[0], positions[0]])
        
        # Case 2: Unbekannte Position
        elif len(positions) == 1 and positions[0] == '?':
            result_rows.extend(generate_combinations(motif, methyl_types[0], enzyme, strands[0]))
        
        # Case 3: Kombination von bekannter und unbekannter Position
        elif len(positions) == 2 and '?' in positions and '-' in strands:
            unknown_index = positions.index('?')
            known_index = 1 - unknown_index
            known_position = positions[known_index]
            known_methyl_type = methyl_types[known_index]
            known_strand = strands[known_index]
            unknown_methyl_type = methyl_types[unknown_index]
            unknown_strand = strands[unknown_index] 

            # Rekursive Analyse der unbekannten Position
            motif_copy = list(motif)
            for i, nucleotide in enumerate(motif_copy):
                if nucleotide == 'C':
                    enzyme_unique_name = f"{enzyme}_C{i+1}"
                    new_motif = motif_copy[:]
                    new_motif[i] = 'X'
                    new_motif = ''.join(new_motif)

                    # Hinzufügen der Kombination für die bekannte Position
                    result_rows.append([row['Enzyme'], enzyme_unique_name, row['Motif'], unknown_methyl_type, unknown_strand, str(i+1)])
                    if known_position.isdigit():
                        result_rows.append([row['Enzyme'], enzyme_unique_name, row['Motif'], known_methyl_type, known_strand, known_position])

        # Case 4: Zwei unbekannte Positionen mit identischem Methyltyp und Strand
        elif len(positions) == 2 and all(pos == '?' for pos in positions) and methyl_types[0] == methyl_types[1] and strands[0] == strands[1]:
            result_rows.extend(generate_combinations(motif, methyl_types[0], enzyme, strands[0]))

        # Case 5: Eine bekannte und eine unbekannte Position auf demselben Strand mit identischem Methyltyp
        elif len(positions) == 2 and '?' in positions and len(set(strands)) == 1:
            unknown_index = positions.index('?')
            known_index = 1 - unknown_index
            
            #Kopie erzeugen für rekursive Analyse
            motif_copy = list(motif)
            
            #Markierung des bekannten Cs mit X damit es aus der rekursiven Analyse fällt
            if positions[known_index].isdigit():
                motif_copy[int(positions[known_index]) - 1] = 'X'
            
            new_motif = ''.join(motif_copy)
            
            #1.Rekursive Suche für restliche Positionen
            new_combinations = generate_combinations(new_motif, methyl_types[unknown_index], enzyme, strands[unknown_index])
            unique_names_set = set()  #duplikate verhindern

            for combination in new_combinations:
                combination[2] = motif  #Motif soll unverändert bleiben
                result_rows.append(combination)
                unique_names_set.add(combination[1])

            # 2.Einträge f+r UniqueEnzymActiity für die bekannten Positonen
            for unique_name in unique_names_set:
                result_rows.append([enzyme, unique_name, motif, methyl_types[known_index], strands[known_index], positions[known_index]])
        
        #Case6 zwei bekannte digit positonen
        elif len(positions) == 2 and positions[0].isdigit() and positions[1].isdigit():
            result_rows.append([enzyme, enzyme, motif, methyl_types[0], strands[0], positions[0]])
            result_rows.append([enzyme, enzyme, motif, methyl_types[1], strands[1], positions[1]])

    # Erstellung eines neuen DataFrames mit den Ergebnissen
    result_df = pd.DataFrame(result_rows, columns=['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif', 'MethylationType', 'Strand', 'Position'])
    return result_df
iupac_to_regex_1 = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'R': '[RAG]',
    'Y': '[YCT]',
    'S': '[SCG]',
    'W': '[WAT]',
    'K': '[KGT]',
    'M': '[MAC]',
    'B': '[BCGTYSK]',
    'D': '[DAGTRWK]',
    'H': '[HACTYWM]',
    'V': '[VACGRSM]',
    'N': '[NATCGRYSWKMBDHVN]'
}

iupac_to_regex_2 = {
    'A': '[AWRMDHVN]',
    'T': '[TYWKBDHVN]',
    'C': '[CYSABDHVN]',
    'G': '[GRSKBDNV]',
    'R': '[AGRSWKMBDHVN]',
    'Y': '[CTYSWKMBDHVN]',
    'S': '[CGSKMBDHVN]',
    'W': '[WATRYKMBDHVN]',
    'K': '[GTKRYSWKBDHVN]',
    'M': '[ACMRYSWMBDHVN]',
    'B': '[CGTRYSWKMBDHVN]',
    'D': '[AGTRYSWKMBDHVN]',
    'H': '[ACTRYSWKMBDHVN]',
    'V': '[ACGRYSWKMBDHVN]',
    'N': '[ATCGRYSWKMBDHVN]'
}

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
    directory= "/Users/azlannisar/Snakemake_Test"
    # Load 5mC Files
    matched_files_5mC = glob.glob(os.path.join(directory, '*_detailed_*_6mA.csv'))
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
    mtase_presence_path = "/Users/azlannisar/Mtase_presence_e_25_values.csv"
    Mtase_presence = pd.read_csv(mtase_presence_path, sep=';', index_col=False)
    
    tsv_file_path = "/Users/azlannisar/TSV_REBASE_data.tsv"
    tsv_df = pd.read_csv(tsv_file_path, sep='\t')

    tsv_df['Enzyme_clean'] = tsv_df['Enzyme'].str.replace('_reverse', '')
    
    tsv_df['NewPosition'] = None
    tsv_df['MethylType'] = None
    tsv_df['Strand'] = None

    # Informationen extrahieren und in neue Spalten speichern
    for i, row in tsv_df.iterrows():
        positions, methyl_types, strands = extract_info(row['Position'])
        if positions[0] is not np.nan:
            tsv_df.at[i, 'NewPosition'] = ','.join(positions)
        else:
            tsv_df.at[i, 'NewPosition'] = np.nan
        if methyl_types[0] is not np.nan:
            tsv_df.at[i, 'MethylType'] = ','.join(methyl_types)
        else:
            tsv_df.at[i, 'MethylType'] = np.nan
        if strands[0] is not np.nan:
            tsv_df.at[i, 'Strand'] = ','.join(strands)
        else:
            tsv_df.at[i, 'Strand'] = np.nan
    #handle the case where we have double modification on the same strand     
    tsv_df.loc[tsv_df['NewPosition'] == '?,?', 'NewPosition'] = '?'
    
    result_df = process_data(tsv_df)
    
    
    
    x_values = [os.path.basename(file_path).split('_detailed_')[1].split('_')[0] for file_path in matched_files_5mC]
    
    enzyme_presence_df = pd.DataFrame(columns=['Enzyme', 'Motif', 'Position', 'Mod_Motif'] + x_values)
    
    
    
    # Populate the DataFrame with e-values and add Mod_Position to Motif column
    for enzyme in Mtase_presence.columns[1:]:  # Skip the 'Isolate' column
        row = [enzyme, np.nan, np.nan, np.nan]  # Initialize Motif and Position columns with NaN
        for x_label in x_values:
            isolate_df = Mtase_presence[Mtase_presence['Isolates'] == x_label]
            if not isolate_df.empty:
                row.append(isolate_df[enzyme].values[0])
            else:
                row.append('NA')
        enzyme_presence_df.loc[len(enzyme_presence_df)] = row
    
    merged_df = pd.merge(enzyme_presence_df, tsv_df[['Enzyme', 'Enzyme_clean', 'Motif', 'Position', 'Mod_Motif']], left_on='Enzyme', right_on='Enzyme_clean', how='left')
    merged_df.rename(columns={'Enzyme_x':'Enzyme', 'Motif_y': 'Motif', 'Position_y': 'Position', 'Mod_Motif_y': 'Mod_Motif'}, inplace=True)
    enzyme_presence_df = merged_df[['Enzyme', 'Motif', 'Position', 'Mod_Motif'] + x_values]
    enzyme_presence_output_path = os.path.join(args.output_dir, "enzyme_presence_df.csv")
    enzyme_presence_df.to_csv(enzyme_presence_output_path, index=False)
    
    #Get unique motifs from the first DataFrame
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
        filtered_positions = enzyme_presence_df[enzyme_presence_df['Position'].str.contains(r'\(6\)')]
        # Plotting one boxplot for each unique motif across all DataFrames
        for motif, scores in motif_scores.items():
            if 'A' in motif:
                clean_motif = motif.replace('"', '')
                regex_pattern_1 = ''.join(iupac_to_regex_1.get(base, base) for base in clean_motif)
                regex_pattern_2 = ''.join(iupac_to_regex_2.get(base, base) for base in clean_motif)
        
                # Attempt to match using the first dictionary (iupac_to_regex_1)
                matching_rows = filtered_positions[filtered_positions['Mod_Motif'].str.contains(regex_pattern_1, na=False)]
        
                # If no match found with first dictionary, try the second (iupac_to_regex_2)
                if matching_rows.empty:
                    matching_rows = filtered_positions[filtered_positions['Mod_Motif'].str.contains(regex_pattern_2, na=False)]
        
                fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the figure size as needed
                box = ax.boxplot(scores, labels=x_values, showfliers=False)
                ax.set_xlabel('')
                ax.set_ylabel('Score')
                ax.set_ylim(0, 110)  # Set the Y-axis limits from 0 to 100
                ax.set_title(f'6mA Boxplot for {motif}')
                plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability 
        
                ax.axhline(y=25, color='magenta', linestyle='--', linewidth=1)
        
                if not matching_rows.empty:
                    # Create a DataFrame for the enzyme presence with e-values
                    table_data = matching_rows.values
                    col_labels = matching_rows.columns
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
                plt.show()