#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:10:21 2024

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
            entry = entry.strip() #White space handling 
            
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
            result_rows_for_unknown = []  # Liste für Kombinationen mit unbekannter Position

            for i, nucleotide in enumerate(motif_copy):
                enzyme_unique_name = None  # Default-Initialisierung
                if unknown_methyl_type == '5mC' or unknown_methyl_type == '4mC':
                    if nucleotide == 'C':
                        enzyme_unique_name = f"{enzyme}_C{i+1}"
                elif unknown_methyl_type == '6mA':
                    if motif.count('A') > 1 and nucleotide == 'A':
                        enzyme_unique_name = f"{enzyme}_A{i+1}"
                    else:
                        enzyme_unique_name = f"{enzyme}"

                if enzyme_unique_name:
                    new_motif = motif_copy[:]
                    new_motif[i] = 'X'
                    new_motif = ''.join(new_motif)

                    # Hinzufügen der Kombination für die unbekannte Position
                    result_rows_for_unknown.append([row['Enzyme'], enzyme_unique_name, row['Motif'], unknown_methyl_type, unknown_strand, str(i+1)])

                    # Hinzufügen der Kombination für die bekannte Position, falls vorhanden
                    if known_position.isdigit():
                        result_rows.append([row['Enzyme'], enzyme_unique_name, row['Motif'], known_methyl_type, known_strand, known_position])

            # Ergebnisse zu result_rows hinzufügen
            result_rows.extend(result_rows_for_unknown)

        # Case 4: Zwei unbekannte Positionen mit identischem Methyltyp und Strand
        elif len(positions) == 2 and all(pos == '?' for pos in positions) and methyl_types[0] == methyl_types[1] and strands[0] == strands[1]:
            result_rows.extend(generate_combinations(motif, methyl_types[0], enzyme, strands[0]))

        # Case 5: Eine bekannte und eine unbekannte Position auf demselben Strand mit identischem Methyltyp
        elif len(positions) == 2 and '?' in positions and len(set(strands)) == 1:
            unknown_index = positions.index('?')
            known_index = 1 - unknown_index
            
            # Kopie erzeugen für rekursive Analyse
            motif_copy = list(motif)
            
            # Markierung des bekannten Cs mit X damit es aus der rekursiven Analyse fällt
            if positions[known_index].isdigit():
                motif_copy[int(positions[known_index]) - 1] = 'X'
            
            new_motif = ''.join(motif_copy)
            
            # 1. Rekursive Suche für restliche Positionen
            new_combinations = generate_combinations(new_motif, methyl_types[unknown_index], enzyme, strands[unknown_index])
            unique_names_set = set()  # Duplikate verhindern

            for combination in new_combinations:
                combination[2] = motif  # Motif soll unverändert bleiben
                result_rows.append(combination)
                unique_names_set.add(combination[1])

            # 2. Einträge für UniqueEnzymActivity für die bekannten Positionen hinzufügen
            for unique_name in unique_names_set:
                result_rows.append([enzyme, unique_name, motif, methyl_types[known_index], strands[known_index], positions[known_index]])
        
        # Case 6: Zwei bekannte Positionen
        elif len(positions) == 2 and positions[0].isdigit() and positions[1].isdigit():
            result_rows.append([enzyme, enzyme, motif, methyl_types[0], strands[0], positions[0]])
            result_rows.append([enzyme, enzyme, motif, methyl_types[1], strands[1], positions[1]])

    # Erstellung eines neuen DataFrames mit den Ergebnissen
    result_df = pd.DataFrame(result_rows, columns=['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif', 'MethylationType', 'Strand', 'Position'])
    return result_df
def read_fasta(file_path):
    contigs = {}
    for record in SeqIO.parse(file_path, "fasta"):
        contigs[record.id] = str(record.seq)
    return contigs
def expand_iupac(motif):
    # Remove any markers (e.g., "")
    motif_cleaned = motif.replace('"', '')
    
    # Expand using iupac_to_regex_2
    expanded_motif = ''.join(iupac_to_regex_2.get(base, base) for base in motif_cleaned)
    
    return expanded_motif


def create_boxplots_with_data(dataframes_list, x_values, contig_positions_dict, enzyme_count_df):
    unique_motifs = dataframes_list[0]['Motif'].unique()

    for motif in unique_motifs:
        if '"C"' not in motif:
            continue

        scores = []
        counts = []
        for df in dataframes_list:
            df_motif = df[df['Motif'] == motif]
            score_values = df_motif['Score1'].values
            scores.append(score_values)
            counts.append(len(score_values))

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.boxplot(scores, labels=x_values, showfliers=False)
        ax.set_xlabel('Isolates')
        ax.set_ylabel('Score')
        ax.set_ylim(0, 110)
        ax.set_title(f'5mC Boxplot for {motif}')
        plt.xticks(rotation=45, ha='right')
        ax.axhline(y=25, color='magenta', linestyle='--', linewidth=1)

        for i, count in enumerate(counts):
            ax.text(i + 1, 10, f'n={count}', ha='center', va='top', fontsize=10, color='black', fontweight='bold')

        collected_data = []
        for isolate in x_values:
            df_isolate = dataframes_list[x_values.index(isolate)]
            motif_data = df_isolate[df_isolate['Motif'] == motif]

            for mod_pos in motif_data['Mod_Pos']:
                if isolate in contig_positions_dict and mod_pos in contig_positions_dict[isolate]:
                    entries = contig_positions_dict[isolate][mod_pos]
                    for entry in entries:
                        if entry['MethylationType'] == '5mC':
                            entry_with_isolate = entry.copy()
                            entry_with_isolate['Isolate'] = isolate
                            collected_data.append(entry_with_isolate)

        if collected_data:
            unique_data = pd.DataFrame(collected_data).drop_duplicates().reset_index(drop=True)

            position_column_index = unique_data.columns.get_loc('Position')
            columns_to_process = unique_data.columns[position_column_index + 1:]

            def update_collected_row(row):
                enzyme = row['Enzyme']
                for identified_column in columns_to_process:
                    if pd.notna(row[identified_column]):  
                        filtered_enzyme_count = enzyme_count_df[enzyme_count_df['Isolate'] == identified_column]
                        filtered_row = filtered_enzyme_count[
                            (filtered_enzyme_count['Enzyme'] == enzyme) &
                            (filtered_enzyme_count['Motif'] == motif)
                        ]

                        if len(filtered_row) > 0:
                            count_value = filtered_row.iloc[0]['Count']
                            row[identified_column] = f"{row[identified_column]} ({count_value})"

                return row

            unique_data = unique_data.apply(update_collected_row, axis=1)
            unique_data = unique_data.loc[:, ~unique_data.columns.isin(['Isolate'])]
            unique_data = unique_data.drop_duplicates().reset_index(drop=True)

            table_data = unique_data.values
            col_labels = unique_data.columns
            cell_text = []

            for row in table_data:
                formatted_row = ['{:.10g}'.format(value) if isinstance(value, (float, np.float64)) else str(value) for value in row]
                cell_text.append(formatted_row)

            table = plt.table(cellText=cell_text, colLabels=col_labels, cellLoc='center', loc='bottom', bbox=[-0.1, -1.7, 1.1, 1.1])
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
            plt.subplots_adjust(left=0.2, bottom=0.6)
        else:
            plt.subplots_adjust(left=0.2, bottom=0.2)

        plt.show()


def count_enzyme_positions_for_motifs_fixed(contig_positions_dict, dataframes_list, x_values):
    enzyme_count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    unique_motifs = dataframes_list[0]['Motif'].unique()
    for motif in unique_motifs:
        for isolate in x_values:
            sp_dict = contig_positions_dict[isolate]
            sp_index = x_values.index(isolate)
            df_sp = dataframes_list[sp_index]

            motif_data = df_sp[df_sp['Motif'] == motif]
            unique_mod_pos = motif_data['Mod_Pos'].unique()

            for mod_pos in unique_mod_pos:
                if mod_pos in sp_dict:
                    entries = sp_dict[mod_pos]
                    for entry in entries:
                        if entry['MethylationType'] == '4mC':
                            enzyme = entry['Enzyme']
                            enzyme_count_dict[isolate][enzyme][motif] += 1

    data = []
    for isolate, enzymes in enzyme_count_dict.items():
        for enzyme, motifs in enzymes.items():
            for motif, count in motifs.items():
                data.append([isolate, enzyme, count, motif])

    count_df = pd.DataFrame(data, columns=['Isolate', 'Enzyme', 'Count', 'Motif'])
    return count_df

iupac_to_regex_2 = {
    'A': '[AWRMDHVN]',
    'T': '[TYWKBDHVN]',
    'C': '[CYSBDHVN]',
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
import re
from collections import defaultdict
from Bio import SeqIO



# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Mtase_presence_Plot1')
    parser.add_argument('MTASE_FILE', type=str, help='Path_to_Mtase_file_path')
    parser.add_argument('Input_dir', type=str, help='Path to list dir')
    parser.add_argument('output_dir', type=str, help='Output directory path')
    parser.add_argument('csv_list', type=str, help='Path to .csv file')
    return parser.parse_args()

def main():
    args = parse_arguments()
    directory = args.Input_dir
    directory= "/Users/azlannisar/Snakemake_Test"
    csv_file_path= args.csv_list
    csv_file_path = '/Users/azlannisar/Data_Test.csv'
    #csv_file_path = "/Users/azlannisar/Snakemake_Test/Data_Test2_updated.csv"
    data_df = pd.read_csv(csv_file_path, names=["File_name", "Reference_path", "pod5_path"], skiprows=1, index_col=False)
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
    mtase_presence_path = "/Users/azlannisar/Mtase_presence_e_25_values.csv"
    Mtase_presence = pd.read_csv(mtase_presence_path, sep=';', index_col=False)
    
    tsv_file_path = "/Users/azlannisar/TSV_REBASE_data.tsv"
    tsv_df = pd.read_csv(tsv_file_path, sep='\t')
    nan_count= tsv_df['Position'].isna().sum()
    print(f"Number of NaN entries in 'Position': {nan_count}")

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
    
    merged_df = pd.merge(enzyme_presence_df, result_df[['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif', 'MethylationType','Strand','Position']], left_on='Enzyme', right_on='Enzyme', how='left')
    merged_df = merged_df.drop(columns=['Motif_x', 'Position_x', 'Mod_Motif'])
    first_columns = ['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif_y', 'MethylationType', 'Strand', 'Position_y']
    other_columns = [col for col in merged_df.columns if col not in first_columns]
    final_columns = first_columns + other_columns
    merged_df = merged_df[final_columns]
    enzyme_presence_df = merged_df.rename(columns={'Motif_y': 'Motif', 'Position_y': 'Position'})
    
    
    enzyme_presence_output_path=os.path.join('/Users/azlannisar/PP_Plot', 'enzyme_presence_df.csv')
    
    enzyme_presence_output_path = os.path.join(args.output_dir, "enzyme_presence_df.csv")
    enzyme_presence_df.to_csv(enzyme_presence_output_path, index=False)
    
    
    
    
    #Create Contig Arrays for each Position we get the Enzymes which are covering the position 
    #Later load this in from arg input csv file row[Reference]

    fasta_files = {}
    for index, row in data_df.iterrows():
        fasta_files[row["File_name"]] = row["Reference_path"]
    
    fasta_files = {
    'SP2565': '/Users/azlannisar/fasta_data/complete_genome/SP2565.fasta',
    'VO31120': '/Users/azlannisar/fasta_data/complete_genome/VO31120.fasta',
    'SP10291': '/Users/azlannisar/fasta_data/complete_genome/SP10291_two_contigs.fasta'
}
    
    #Read Fasta Files and create contig_positions dict precreate the array size 
    contig_positions_dict = defaultdict(lambda: defaultdict(list))

    # Iterate over each isolate
    for isolate in x_values:
        contig_sequences = read_fasta(fasta_files[isolate])
        
        # Keep track of positions across all contigs for the current isolate
        isolate_positions = defaultdict(list)
        
        # Iterate over each contig sequence in the current isolate's FASTA file
        for contig_name, contig_sequence in contig_sequences.items():
            # Iterate through enzyme_presence dataframe
            for rowI in range(len(enzyme_presence_df)):
                # Extract complete row as data
                row_data = enzyme_presence_df.loc[rowI].to_dict()
                
                motif = row_data['Motif']
                position = int(row_data['Position'])
                strand = row_data.get('Strand', '+')  # Assuming there's a 'Strand' column in your dataframe
                
                # Transformation
                cleaned_motif = motif.replace('"', '')  # Remove markers
                expanded_motif = expand_iupac(cleaned_motif)
                
                # Compile the regular expression as motif
                p = re.compile(expanded_motif)
                
                # Find all matches of the regular expression motif in contig_sequence
                for m in p.finditer(contig_sequence):
                    match_start_pos = m.start()
                    
                    if strand == "+":
                        match_methylation_pos = match_start_pos + position - 1
                    elif strand == "-":
                        # Calculate contig position from the right for reverse strand
                        match_methylation_pos = match_start_pos + len(motif) - position
                    
                    if 0 <= match_methylation_pos < len(contig_sequence):
                        # Append data to the matching isolate and contig position if not already present
                        if row_data not in isolate_positions[match_methylation_pos]:
                            isolate_positions[match_methylation_pos].append(row_data)
        
        # Merge isolate_positions into the main contig_positions_dict for the isolate
        for pos, data_list in isolate_positions.items():
            for data in data_list:
                if data not in contig_positions_dict[isolate][pos]:
                    contig_positions_dict[isolate][pos].append(data)
    
    enzyme_position_count_df = count_enzyme_positions_for_motifs_fixed(contig_positions_dict, dataframes_list, x_values)
    
    #Get unique motifs from the first DataFrame
    create_boxplots_with_data(dataframes_list, x_values, contig_positions_dict, enzyme_position_count_df)