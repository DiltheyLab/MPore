#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 15:53:18 2024

@author: azlannisar
"""

import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from Bio import SeqIO
import argparse

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Mtase_presence_Plot1')
    parser.add_argument('MTASE_FILE', type=str, help='Path to Mtase file')
    parser.add_argument('Input_dir', type=str, help='Path to input directory')
    parser.add_argument('output_dir', type=str, help='Output directory path')
    parser.add_argument('csv_list', type=str, help='Path to .csv file')
    return parser.parse_args()

def read_fasta(file_path):
    contigs = {}
    for record in SeqIO.parse(file_path, "fasta"):
        contigs[record.id] = str(record.seq)
    return contigs

def extract_info(position):
    if pd.isna(position):
        return [np.nan], [np.nan], [np.nan]
    else:
        position = str(position)
        positions = []
        methyl_types = []
        strands = []
        entries = position.split(',')
        for entry in entries:
            entry = entry.strip()
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

        if len(positions) == 1 and positions[0] != '?':
            result_rows.append([enzyme, enzyme, motif, methyl_types[0], strands[0], positions[0]])
        elif len(positions) == 1 and positions[0] == '?':
            result_rows.extend(generate_combinations(motif, methyl_types[0], enzyme, strands[0]))
        elif len(positions) == 2 and '?' in positions and '-' in strands:
            unknown_index = positions.index('?')
            known_index = 1 - unknown_index
            known_position = positions[known_index]
            known_methyl_type = methyl_types[known_index]
            known_strand = strands[known_index]
            unknown_methyl_type = methyl_types[unknown_index]
            unknown_strand = strands[unknown_index]

            motif_copy = list(motif)
            result_rows_for_unknown = []

            for i, nucleotide in enumerate(motif_copy):
                enzyme_unique_name = None
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
                    result_rows_for_unknown.append([row['Enzyme'], enzyme_unique_name, row['Motif'], unknown_methyl_type, unknown_strand, str(i+1)])
                    if known_position.isdigit():
                        result_rows.append([row['Enzyme'], enzyme_unique_name, row['Motif'], known_methyl_type, known_strand, known_position])

            result_rows.extend(result_rows_for_unknown)

        elif len(positions) == 2 and all(pos == '?' for pos in positions) and methyl_types[0] == methyl_types[1] and strands[0] == strands[1]:
            result_rows.extend(generate_combinations(motif, methyl_types[0], enzyme, strands[0]))

        elif len(positions) == 2 and '?' in positions and len(set(strands)) == 1:
            unknown_index = positions.index('?')
            known_index = 1 - unknown_index
            motif_copy = list(motif)
            if positions[known_index].isdigit():
                motif_copy[int(positions[known_index]) - 1] = 'X'
            new_motif = ''.join(motif_copy)
            new_combinations = generate_combinations(new_motif, methyl_types[unknown_index], enzyme, strands[unknown_index])
            unique_names_set = set()

            for combination in new_combinations:
                combination[2] = motif
                result_rows.append(combination)
                unique_names_set.add(combination[1])

            for unique_name in unique_names_set:
                result_rows.append([enzyme, unique_name, motif, methyl_types[known_index], strands[known_index], positions[known_index]])

        elif len(positions) == 2 and positions[0].isdigit() and positions[1].isdigit():
            result_rows.append([enzyme, enzyme, motif, methyl_types[0], strands[0], positions[0]])
            result_rows.append([enzyme, enzyme, motif, methyl_types[1], strands[1], positions[1]])

    result_df = pd.DataFrame(result_rows, columns=['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif', 'MethylationType', 'Strand', 'Position'])
    return result_df

def expand_iupac(motif):
    motif_cleaned = motif.replace('"', '')
    expanded_motif = ''.join(iupac_to_regex_2.get(base, base) for base in motif_cleaned)
    return expanded_motif

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
                        if entry['MethylationType'] == '6mA':
                            enzyme = entry['Enzyme']
                            enzyme_count_dict[isolate][enzyme][motif] += 1

    data = []
    for isolate, enzymes in enzyme_count_dict.items():
        for enzyme, motifs in enzymes.items():
            for motif, count in motifs.items():
                data.append([isolate, enzyme, count, motif])

    count_df = pd.DataFrame(data, columns=['Isolate', 'Enzyme', 'Count', 'Motif'])
    return count_df

def create_boxplots_with_data(dataframes_list, x_values, contig_positions_dict, enzyme_count_df, output_dir):
    unique_motifs = dataframes_list[0]['Motif'].unique()

    for motif in unique_motifs:
        if '"A"' not in motif:
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
        ax.set_title(f'6mA Boxplot for {motif}')
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
                        if entry['MethylationType'] == '6mA':
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
                            row[identified_column] = f"{row[identified_column]:.2e} ({count_value})"

                return row

            unique_data = unique_data.apply(update_collected_row, axis=1)
            unique_data = unique_data.loc[:, ~unique_data.columns.isin(['Isolate'])]
            unique_data = unique_data.drop_duplicates().reset_index(drop=True)

            table_data = unique_data.values
            col_labels = unique_data.columns
            cell_text = []

            for row in table_data:
                formatted_row = ['{:.2e}'.format(value) if isinstance(value, (float, np.float64)) else str(value) for value in row]
                cell_text.append(formatted_row)

            table = plt.table(cellText=cell_text, colLabels=col_labels, cellLoc='center', loc='bottom', bbox=[-0.1, -1.7, 1.1, 1.1])
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.5)
            plt.subplots_adjust(left=0.2, bottom=0.6)
        else:
            plt.subplots_adjust(left=0.2, bottom=0.2)

        plt.savefig(os.path.join(output_dir, f'{motif}_6mA_boxplot.png'))
        plt.close()

        


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

def main():
    args = parse_arguments()
    directory = args.Input_dir
    directory = "/Users/azlannisar/Snakemake_Output_final/Output2"
    csv_file_path= "/Users/azlannisar/Snakemake_Output_final/Data_Test.csv"
    csv_file_path = args.csv_list

    # Load CSV list
    data_df = pd.read_csv(csv_file_path, names=["File_name", "Reference_path", "pod5_path"], skiprows=1, index_col=False)

    # Load 5mC Files
    matched_files_5mC = glob.glob(os.path.join(directory, '*_detailed_*_6mA.csv'))
    dataframes_list = []
    for file_path in matched_files_5mC:
        try:
            df = pd.read_csv(file_path)
            df = df[df['Score1'] != '/']
            df['Score1'] = pd.to_numeric(df['Score1'], errors='coerce')
            dataframes_list.append(df)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            
    
    

    # Read Mtase_presence DataFrame
    mtase_presence_path = args.MTASE_FILE
    mtase_presence_path = "/Users/azlannisar/Snakemake_Output_final/Output2/Mtase_presence_e_25_values.csv"
    Mtase_presence = pd.read_csv(mtase_presence_path, sep=',', index_col=False)

    # Read TSV file
    tsv_file_path = "/Users/azlannisar/TSV_REBASE_data.tsv"
    tsv_df = pd.read_csv(tsv_file_path, sep='\t')
    nan_count = tsv_df['Position'].isna().sum()
    print(f"Number of NaN entries in 'Position': {nan_count}")

    tsv_df['Enzyme_clean'] = tsv_df['Enzyme'].str.replace('_reverse', '')
    tsv_df['NewPosition'] = None
    tsv_df['MethylType'] = None
    tsv_df['Strand'] = None

    for i, row in tsv_df.iterrows():
        positions, methyl_types, strands = extract_info(row['Position'])
        tsv_df.at[i, 'NewPosition'] = ','.join(positions) if positions[0] is not np.nan else np.nan
        tsv_df.at[i, 'MethylType'] = ','.join(methyl_types) if methyl_types[0] is not np.nan else np.nan
        tsv_df.at[i, 'Strand'] = ','.join(strands) if strands[0] is not np.nan else np.nan

    tsv_df.loc[tsv_df['NewPosition'] == '?,?', 'NewPosition'] = '?'
    result_df = process_data(tsv_df)

    x_values = [os.path.basename(file_path).split('_detailed_')[1].split('_')[0] for file_path in matched_files_5mC]

    enzyme_presence_df = pd.DataFrame(columns=['Enzyme', 'Motif', 'Position', 'Mod_Motif'] + x_values)

    for enzyme in Mtase_presence.columns[1:]:
        row = [enzyme, np.nan, np.nan, np.nan]
        for x_label in x_values:
            isolate_df = Mtase_presence[Mtase_presence['Isolates'] == x_label]
            row.append(isolate_df[enzyme].values[0] if not isolate_df.empty else 'NA')
        enzyme_presence_df.loc[len(enzyme_presence_df)] = row

    merged_df = pd.merge(enzyme_presence_df, result_df[['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif', 'MethylationType', 'Strand', 'Position']], on='Enzyme', how='left')
    merged_df = merged_df.drop(columns=['Motif_x', 'Position_x', 'Mod_Motif'])
    first_columns = ['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif_y', 'MethylationType', 'Strand', 'Position_y']
    other_columns = [col for col in merged_df.columns if col not in first_columns]
    final_columns = first_columns + other_columns
    merged_df = merged_df[final_columns]
    enzyme_presence_df = merged_df.rename(columns={'Motif_y': 'Motif', 'Position_y': 'Position'})

    enzyme_presence_output_path = os.path.join(args.output_dir, "enzyme_presence_df.csv")
    enzyme_presence_df.to_csv(enzyme_presence_output_path, index=False)

    fasta_files = {}
    for index, row in data_df.iterrows():
        fasta_files[row["File_name"]] = row["Reference_path"]

    contig_positions_dict = defaultdict(lambda: defaultdict(list))

    for isolate in x_values:
        contig_sequences = read_fasta(fasta_files[isolate])
        isolate_positions = defaultdict(list)

        for contig_name, contig_sequence in contig_sequences.items():
            for rowI in range(len(enzyme_presence_df)):
                row_data = enzyme_presence_df.loc[rowI].to_dict()
                motif = row_data['Motif']
                position = int(row_data['Position'])
                strand = row_data.get('Strand', '+')
                cleaned_motif = motif.replace('"', '')
                expanded_motif = expand_iupac(cleaned_motif)
                p = re.compile(expanded_motif)

                for m in p.finditer(contig_sequence):
                    match_start_pos = m.start()
                    if strand == "+":
                        match_methylation_pos = match_start_pos + position - 1
                    elif strand == "-":
                        match_methylation_pos = match_start_pos + len(motif) - position

                    if 0 <= match_methylation_pos < len(contig_sequence):
                        if row_data not in isolate_positions[match_methylation_pos]:
                            isolate_positions[match_methylation_pos].append(row_data)

        for pos, data_list in isolate_positions.items():
            for data in data_list:
                if data not in contig_positions_dict[isolate][pos]:
                    contig_positions_dict[isolate][pos].append(data)

    enzyme_position_count_df = count_enzyme_positions_for_motifs_fixed(contig_positions_dict, dataframes_list, x_values)
    
    #Load Bed Files into list for readbases analysis 
    matched_files_bed = glob.glob(os.path.join(directory, '*.bed'))
    dataframes_list_bed = []
    column_names = ['sample', 'start', 'end', 'name', 'score', 'strand', 'tstart', 'tend', 'color', 'coverage', 'frequence', 'modified', 'canonical', 'other_mod','delted','failed','substitution','No_call']
    for file_path in matched_files_bed:
        try:
            df = pd.read_csv(file_path, sep= '\t', header=None, names=column_names)
            df = df[df['name']== 'a']
            df['score'] = pd.to_numeric(df['score'], errors='coerce')
            dataframes_list_bed.append(df)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    
    
    #Create Big_List with the data 
    Big_list = []
    for df in dataframes_list:
    
        new_df = pd.DataFrame({
            'Base': 'A', 
            'pos': df['Mod_Pos'],
            'strand':df['Strand']
            
        })
        
        # Append the new dataframe to Big_list
        Big_list.append(new_df)
    #Create Dataframes with Reads based and Pos Based data
    updated_big_list = []
    assert len(dataframes_list_bed) == len(Big_list), "The lists must have the same number of elements."
    for df_bed, df_big in zip(dataframes_list_bed, Big_list):
        try:
            df_bed['end'] = pd.to_numeric(df_bed['end'], errors='coerce')
            df_big['pos'] = pd.to_numeric(df_big['pos'], errors='coerce')
    
            merged_df = pd.merge(df_big, df_bed[['end', 'coverage', 'modified']],
                                 left_on='pos', right_on='end', how='left')
            
            
            merged_df.drop(columns=['end'], inplace=True)
            merged_df.dropna(subset=['coverage', 'modified'], inplace=True)
    
            # Append the updated dataframe to the new list
            updated_big_list.append(merged_df)
        except Exception as e:
            print(f"Error processing dataframes: {e}")
    updated_dict = {}
    assert len(x_values) == len(updated_big_list), "x_values and updated_big_list must have the same number of elements."
    for key, df in zip(x_values, updated_big_list):
        updated_dict[key] = df
    
    #Add the Indikators for each Dataframe
    for key, updated_df in updated_dict.items():
        contig_data = contig_positions_dict.get(key, defaultdict(list))
        unique_enzymes = set()
        
        #Get unique enzyme_names from contig_positions_dict
        for enzyme_info_list in contig_data.values():
            for enzyme_info in enzyme_info_list:
                enzyme_name = enzyme_info.get('Enzyme')
                if enzyme_name:
                    unique_enzymes.add(enzyme_name)
        
        #Create columns based on unique_names and 0 initialization 
        for enzyme in unique_enzymes:
            updated_df[enzyme] = 0
            
        for Position, enzyme_info_list in contig_data.items():
            for enzyme_info in enzyme_info_list:
                enzyme_name = enzyme_info.get('Enzyme')
                if enzyme_name and enzyme_name in updated_df.columns:
                    matching_indices = updated_df[updated_df['pos'] == Position].index
                    updated_df.loc[matching_indices, enzyme_name] = 1
        updated_dict[key] = updated_df
        
    output_dir = "/Users/azlannisar/Snakemake_Output_final"
    os.makedirs(output_dir, exist_ok=True)
    for key, df in updated_dict.items():
        filename= os.path.join(output_dir, f"{key}.csv")
        df.to_csv(filename, index = False)

    #Save as complete dict .json
    import json
    serializable_dict = {key: df.to_dict(orient='records') for key, df in updated_dict.items()}
    with open(os.path.join(output_dir, 'updated_dict.json'), 'w') as f:
        json.dump(serializable_dict, f)

    

    #output_test = args.output_dir
    #output_test= "/Users/azlannisar/Snakemake_Output_final/Output/plots_4mC"
    create_boxplots_with_data(dataframes_list, x_values, contig_positions_dict, enzyme_position_count_df, args.output_dir)

if __name__ == "__main__":
    main()