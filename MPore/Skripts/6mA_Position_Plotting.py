#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 18:14:31 2024

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
from Bio.Seq import Seq
import argparse



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
            complement_motif = str(Seq(motif).complement())[::-1]
            if complement_motif != motif:
                complement_results = generate_combinations(complement_motif, methyl_types[0], enzyme, '-')
                for result in complement_results:
                    result[2] = motif
                result_rows.extend(complement_results)
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
    expanded_motif = ''.join(iupac_to_regex.get(base, base) for base in motif_cleaned)
    return expanded_motif

def count_enzyme_positions_for_motifs_fixed(contig_positions_dict, dataframes_list, x_values):
    
    enzyme_count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    unique_motifs = dataframes_list[0]['Motif'].unique()
        
    for motif in unique_motifs:
        
        for isolate in x_values:
            
            sp_dict = contig_positions_dict[isolate]
            sp_index = x_values.index(isolate)
            df_sp = dataframes_list[sp_index]

            
            motif_data = df_sp[df_sp['Motif'] == motif].copy()
            motif_data['Contig'] = motif_data['Contig'].astype(str)
            unique_mod_pos = motif_data[['Mod_Pos', 'Contig']].drop_duplicates()

            for mod_pos, contig in unique_mod_pos.values:
                
                if (str(contig), mod_pos) in sp_dict:
                    entries = sp_dict[(str(contig), mod_pos)]
                                        
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

def create_boxplots_with_data(dataframes_list, x_values, contig_positions_dict, enzyme_count_df, log_results_df, output_dir):
    unique_motifs = dataframes_list[0]['Motif'].unique()
    
    
    for motif in unique_motifs:
        if '"A"' not in motif:
            continue

        scores = []
        counts = []
        for df in dataframes_list:
            df_motif = df[df['Motif'] == motif]
            scores.append(df_motif['Score1'].values)
            counts.append(len(df_motif))

        collected_data = []
        for isolate in x_values:
            df_isolate = dataframes_list[x_values.index(isolate)]
            motif_data = df_isolate[df_isolate['Motif'] == motif]

            for mod_pos, contig in zip(motif_data['Mod_Pos'], motif_data['Contig']):
                key = (str(contig), int(mod_pos))
                if isolate in contig_positions_dict and key in contig_positions_dict[isolate]:
                    for entry in contig_positions_dict[isolate][key]:
                        if entry['MethylationType'] == '6mA':
                            new_entry = entry.copy()
                            new_entry['Isolate'] = isolate
                            new_entry.pop('Contig', None)
                            collected_data.append(new_entry)

        fig_width = max(12, len(x_values) * 1.5)
        num_enzyme_rows = len(enzyme_count_df['Enzyme'].unique())
        fig_height = max(10, num_enzyme_rows * 0.6)

        fig, (ax_box, ax_table) = plt.subplots(
            2, 1,
            figsize=(fig_width, fig_height),
            gridspec_kw={'height_ratios': [2, 1]},
            constrained_layout=True
        )

        # === Boxplot ===
        ax_box.boxplot(scores, labels=x_values, showfliers=False)
        ax_box.set_xlabel('Isolates', fontsize=14)
        ax_box.set_ylabel('Score', fontsize=14)
        ax_box.set_ylim(0, 110)
        ax_box.set_title(f'6mA Boxplot for {motif}', fontsize=16)
        ax_box.axhline(y=25, color='magenta', linestyle='--', linewidth=1)
        ax_box.tick_params(axis='x', rotation=45)

        for i, count in enumerate(counts):
            ax_box.text(i + 1, 10, f'n={count}', ha='center', va='top', fontsize=10)

        # === Table ===
        if collected_data:
            unique_data = pd.DataFrame(collected_data).drop_duplicates().reset_index(drop=True)
            unique_data.rename(columns={'EnzymeUniqueMethylationActivity': 'Uniquemethylation'}, inplace=True)

            pos_idx = unique_data.columns.get_loc('Position')
            columns_to_process = unique_data.columns[pos_idx + 1:]

            def update_row(row):
                enzyme = row['Enzyme']
                for iso in columns_to_process:
                    if pd.notna(row[iso]):
                        count_row = enzyme_count_df[
                            (enzyme_count_df['Isolate'] == iso) &
                            (enzyme_count_df['Enzyme'] == enzyme) &
                            (enzyme_count_df['Motif'] == motif)
                        ]
                        count_val = count_row.iloc[0]['Count'] if not count_row.empty else 'NA'
                        
                        try:
                            val = float(row[iso])
                            formatted_val = f"{val:.2e}"
                        except (ValueError, TypeError):
                            val = row[iso]
                            formatted_val = f"{val}"
                        
                        beta_row = log_results_df[
                            (log_results_df['Isolate'].str.lower() == iso.lower()) &
                            (log_results_df['Enzyme'].str.lower() == enzyme.lower())
                        ]
                        if not beta_row.empty:
                            beta_val = beta_row.iloc[0]['beta_coefficient']
                            row[iso] = f"{formatted_val} ({count_val}), ß={beta_val:.2f}"
                        else:
                            row[iso] = f"{formatted_val}({count_val})"
                return row

            unique_data = unique_data.apply(update_row, axis=1)
            unique_data = unique_data.drop(columns=['Isolate'], errors='ignore')
            unique_data = unique_data.drop_duplicates().reset_index(drop=True)
            table_data = unique_data.values
            col_labels = unique_data.columns

            ax_table.axis('tight')
            ax_table.axis('off')
            table = ax_table.table(
                cellText=table_data,
                colLabels=col_labels,
                cellLoc='center',
                loc='center'
            )

            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(2.0, 1.5)
        else:
            ax_table.text(0.5, 0.5, f'No data for motif: {motif}', ha='center', va='center', fontsize=12)
            ax_table.axis('off')

        filename = os.path.join(output_dir, f'{motif}_6mA_boxplot.png')
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        
iupac_to_regex = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[CG]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
    'N': '[ATCG]'
}

def main():
    args = parse_arguments()
    directory = args.Input_dir
    csv_file_path = args.csv_list
    data_df = pd.read_csv(csv_file_path, names=["File_name", "Reference_path", "pod5_path"], skiprows=1, index_col=False)
    file_names = set(data_df["File_name"].astype(str))

    matched_files_5mC = glob.glob(os.path.join(directory, '*_detailed_*_6mA.csv'))
    dataframes_list = []
    for file_path in matched_files_5mC:
        base_name = os.path.basename(file_path)
        if any(f"Sample_DF_detailed_{file_name}_6mA.csv" == base_name for file_name in file_names):
            try:
                df = pd.read_csv(file_path, low_memory=False)
                df = df[df['Score1'] != '/']
                df['Score1'] = pd.to_numeric(df['Score1'], errors='coerce')
                dataframes_list.append(df)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
    
    
    mtase_presence_path = args.MTASE_FILE
    
    Mtase_presence = pd.read_csv(mtase_presence_path, sep=',', index_col=False)
    Mtase_presence = Mtase_presence[Mtase_presence["Isolates"].isin(file_names)]
    Mtase_presence = Mtase_presence.dropna(axis=1, how='all')

    matched_beta_csv=[]
    matched_beta_csv= glob.glob(os.path.join(directory, 'Beta_coef_p_values_filt_6mA_*.csv'))
    Log_Results_df=pd.DataFrame()
    for beta in matched_beta_csv:
        try:
            beta_df=pd.read_csv(beta)
            Log_Results_df=pd.concat([Log_Results_df, beta_df], ignore_index=True)
        except Exception as e:
            print(f"Error reading {beta}: {e}")
    
    Log_Results_df = Log_Results_df[Log_Results_df['Isolate'].isin(file_names)]
        
    
    tsv_file_path="./TSV_Enzyme.csv"
    
    tsv_df = pd.read_csv(tsv_file_path, sep=',')
    nan_count = tsv_df['Position'].isna().sum()
    print(f"Number of NaN entries in 'Position': {nan_count}")
    tsv_df = tsv_df[['Enzyme', 'Motif', 'Position']].copy()
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

    x_values = [
    file_name
    for file_path in matched_files_5mC
        if (file_name := os.path.basename(file_path).split('_detailed_')[1].split('_6mA')[0]) in file_names
    ]

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
    enzyme_presence_df = enzyme_presence_df.dropna(subset=['Position']).reset_index(drop=True)
    enzyme_presence_df['Motif'] = enzyme_presence_df['Motif'].str.extract(r'([A-Z]+)')

    enzyme_presence_output_path = os.path.join(args.output_dir, "enzyme_presence_df.csv")
    enzyme_presence_df.to_csv(enzyme_presence_output_path, index=False)
    Set_of_non_negatives= Log_Results_df[Log_Results_df['beta_coefficient'] >= 0]
    enzyme_presence_df=enzyme_presence_df[enzyme_presence_df['Enzyme'].isin(Set_of_non_negatives['Enzyme'])]    

    fasta_files = {}
    for index, row in data_df[data_df["File_name"].isin(file_names)].iterrows():
        fasta_files[row["File_name"]] = row["Reference_path"]
            
    
    contig_positions_dict = defaultdict(lambda: defaultdict(list))
    
    
    for isolate in x_values:
        contig_sequences = read_fasta(fasta_files[isolate])  
        isolate_positions = defaultdict(list)

        seen_entries = defaultdict(set)

        for contig_name, contig_sequence in contig_sequences.items():
            
            for rowI in range(len(enzyme_presence_df)):
                
                row_data = enzyme_presence_df.iloc[rowI].to_dict()
                motif = row_data['Motif']
                position = int(row_data['Position'])
                strand = row_data['Strand']
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
                        
                        row_data_with_contig = {**row_data, 'Contig': contig_name}
 
                        entry_key = (
                            row_data_with_contig['EnzymeUniqueMethylationActivity'],
                            row_data_with_contig['Motif'],
                            match_methylation_pos,
                        )
    
                        if entry_key not in seen_entries[contig_name, match_methylation_pos]:
                            isolate_positions[(contig_name, match_methylation_pos)].append(row_data_with_contig)
                            seen_entries[contig_name, match_methylation_pos].add(entry_key)
    
        
        for (contig_name, pos), data_list in isolate_positions.items():
            contig_positions_dict[isolate][(contig_name, pos)].extend(data_list)

    contig_positions_dict_1 = defaultdict(lambda: defaultdict(list))
    
    
    for isolate in x_values:
        contig_sequences = read_fasta(fasta_files[isolate])  
        isolate_positions = defaultdict(list)

        seen_entries = defaultdict(set)

        for contig_name, contig_sequence in contig_sequences.items():
            
            for rowI in range(len(enzyme_presence_df)):
                
                row_data = enzyme_presence_df.iloc[rowI].to_dict()
                motif = row_data['Motif']
                position = int(row_data['Position'])
                strand = row_data['Strand']

                reverse_complement_motif = str(Seq(motif).reverse_complement())

                cleaned_motif = reverse_complement_motif.replace('"', '')
                expanded_motif = expand_iupac(cleaned_motif)  
                p = re.compile(expanded_motif)
    
                
                for m in p.finditer(contig_sequence):
                    match_start_pos = m.start()
                    
                    if strand == "+":
                        match_methylation_pos = match_start_pos + len(motif) - position
                    elif strand == "-":
                        match_methylation_pos = match_start_pos + position -1 
    
                    
                    if 0 <= match_methylation_pos < len(contig_sequence):
                        
                        row_data_with_contig = {**row_data, 'Contig': contig_name}
    
                        
                        entry_key = (
                            row_data_with_contig['EnzymeUniqueMethylationActivity'],
                            match_methylation_pos,
                        )
    
                        
                        if entry_key not in seen_entries[contig_name, match_methylation_pos]:
                            isolate_positions[(contig_name, match_methylation_pos)].append(row_data_with_contig)
                            seen_entries[contig_name, match_methylation_pos].add(entry_key)
    
        
        for (contig_name, pos), data_list in isolate_positions.items():
            contig_positions_dict_1[isolate][(contig_name, pos)].extend(data_list)
        
    combined_contig_positions_dict = defaultdict(lambda: defaultdict(list))

    for isolate, contig_data in contig_positions_dict.items():
        for (contig_name, pos), data_list in contig_data.items():
            
            combined_contig_positions_dict[isolate][(contig_name, pos)].extend(data_list)

    
    for isolate, contig_data in contig_positions_dict_1.items():
        for (contig_name, pos), data_list in contig_data.items():
            
            combined_contig_positions_dict[isolate][(contig_name, pos)].extend(data_list) 
    enzyme_position_count_df = count_enzyme_positions_for_motifs_fixed(combined_contig_positions_dict, dataframes_list, x_values)

    
    
    create_boxplots_with_data(dataframes_list, x_values, combined_contig_positions_dict, enzyme_position_count_df, Log_Results_df, args.output_dir)

if __name__ == "__main__":
    main()


