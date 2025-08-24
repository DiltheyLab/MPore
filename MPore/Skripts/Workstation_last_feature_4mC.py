#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:47:23 2024

@author: azlannisar
"""

import os
import glob
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import json
import math

def is_palindrome(motif):
    return motif == str(Seq(motif).reverse_complement())

def read_fasta(file_path):
    contigs = {}
    for record in SeqIO.parse(file_path, "fasta"):
        contigs[record.id] = str(record.seq)
    return contigs

def expand_iupac(motif):
    motif_cleaned = motif.replace('"', '')
    expanded_motif = ''.join(iupac_to_regex_2.get(base, base) for base in motif_cleaned)
    return expanded_motif

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

def parse_arguments():
    parser = argparse.ArgumentParser(description='Data refinement for log analysis')
    parser.add_argument('MTASE_FILE', type=str, help='Path to Mtase file')
    parser.add_argument('Input_dir', type=str, help='Path to input directory')
    parser.add_argument('output_dir', type=str, help='Output directory path')
    parser.add_argument('csv_list', type=str, help='Path to .csv file')
    return parser.parse_args()
def main():
    args = parse_arguments()
    directory=args.Input_dir
    csv_file_path=args.csv_list
    df = pd.read_csv(csv_file_path)
    file_names = df['File_name'].astype(str).tolist()
    matched_files_5mC = []
    matched_files_5mC = [
        file for file in glob.glob(os.path.join(directory, '*_detailed_*_4mC.csv'))
        if any(f"Sample_DF_detailed_{file_name}_4mC.csv" in file for file_name in file_names)]
    dataframes_list = []
    matched_files_5mC.sort()
    for file_path in matched_files_5mC:
        try:
            df = pd.read_csv(file_path, low_memory=False)
            df = df[df['Score1'] != '/']
            df['Score1'] = pd.to_numeric(df['Score1'], errors='coerce')
            dataframes_list.append(df)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    matched_files_bed = [
        file for file in glob.glob(os.path.join(directory, '*.bed'))
        if any(f"{file_name}.bed" in file for file_name in file_names)
    ]
    dataframes_list_bed = []
    matched_files_bed.sort()
    column_names = ['sample', 'start', 'end', 'name', 'score', 'strand', 'tstart', 'tend', 'color', 'coverage', 'frequence', 'modified', 'canonical', 'other_mod','delted','failed','substitution','No_call']
    for file_path in matched_files_bed:
        try:
            df = pd.read_csv(file_path, sep= '\t', header=None, names=column_names)
            df = df[df['name']== '21839']
            df['score'] = pd.to_numeric(df['score'], errors='coerce')
            dataframes_list_bed.append(df)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    
    mtase_presence_path=args.MTASE_FILE
    Mtase_presence = pd.read_csv(mtase_presence_path, sep=',', index_col=False)
    Mtase_presence = Mtase_presence.dropna(axis=1, how='all')


    
    Big_list = []
    for df in dataframes_list:
    
        new_df = pd.DataFrame({
            'Base': 'C', 
            'pos': df['Mod_Pos'],
            'strand':df['Strand'],
            'score':df['Score1'],
            'motif':df['Motif'],
            'contig':df['Contig']
            
        })
        
        Big_list.append(new_df)
    
    updated_big_list = []
    assert len(dataframes_list_bed) == len(Big_list) 
    for df_bed, df_big in zip(dataframes_list_bed, Big_list):
        try:
            df_bed['start'] = pd.to_numeric(df_bed['start'], errors='coerce')
            df_big['pos'] = pd.to_numeric(df_big['pos'], errors='coerce')
    
            merged_df = pd.merge(df_big, df_bed[['start', 'coverage', 'modified','sample']],
                                 left_on=['pos', 'contig'], right_on=['start','sample'], how='left')
            
            
            merged_df.drop(columns=['start'], inplace=True)
            merged_df.dropna(subset=['coverage', 'modified'], inplace=True)

            updated_big_list.append(merged_df)
        except Exception as e:
            print(f"Error processing dataframes: {e}")
    
    data_df = pd.read_csv(csv_file_path, names=["File_name", "Reference_path", "pod5_path"], skiprows=1, index_col=False)
    data_df = data_df.sort_values(by="File_name").reset_index(drop=True)

    tsv_file_path = "./TSV_Enzyme.csv"
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
    x_values = data_df['File_name'].dropna().tolist()
    
    
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
    merged_df = merged_df.dropna(subset=[merged_df.columns[1]])
    enzyme_presence_df = merged_df.rename(columns={'Motif_y': 'Motif', 'Position_y': 'Position'})
    enzyme_presence_df = enzyme_presence_df.dropna(subset=['Position']).reset_index(drop=True)

    fasta_files = {}
    for index, row in data_df.iterrows():
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

    updated_dict = {}
    assert len(x_values) == len(updated_big_list), "x_values and updated_big_list must have the same number of elements."
    for key, df in zip(x_values, updated_big_list):
        updated_dict[key] = df
    
    for key, updated_df in updated_dict.items():
        combined_contig_data = combined_contig_positions_dict.get(key, defaultdict(list)) 
        unique_enzymes = set()
        unique_enzymes = {
            enzyme_info['Enzyme']
            for enzyme_info_list in combined_contig_data.values()
            for enzyme_info in enzyme_info_list
            if enzyme_info.get('Enzyme')
        }
        
        unique_enzymes = list(unique_enzymes)
        enzyme_columns = pd.DataFrame(0, index=updated_df.index, columns=unique_enzymes)

        
        updated_df = pd.concat([updated_df, enzyme_columns], axis=1)        
        
        lookup_no_strand = defaultdict(set)
        for contig_number, enzyme_info_list in combined_contig_data.items():
            for enzyme_info in enzyme_info_list:
                enzyme_name = enzyme_info.get('Enzyme')
                methylation_type = enzyme_info.get('MethylationType')
                pos = enzyme_info.get('Position')  

        
                if enzyme_name and methylation_type == '4mC' and pos is not None:
                    lookup_no_strand[(contig_number)].add(enzyme_name)
    
        for idx, row in updated_df.iterrows():
            contig = row['contig']
            pos = row['pos']
            enzymes = lookup_no_strand.get((contig, pos), set())
            for enzyme in enzymes:
                updated_df.at[idx, enzyme] = 1
    
        updated_dict[key] = updated_df

        non_enzymatic_cols = ['Base', 'pos', 'strand', 'motif', 'contig', 'coverage', 'modified', 'sample']
    
        for key, updated_df in updated_dict.items():
            enzyme_cols = [col for col in updated_df.columns if col not in non_enzymatic_cols]
            cols_to_drop = [col for col in enzyme_cols if updated_df[col].sum() == 0]
            updated_df.drop(columns=cols_to_drop, inplace=True)
            updated_dict[key] = updated_df


    output_dir=args.output_dir
    serializable_dict = {key: df.to_dict(orient='records') for key, df in updated_dict.items()}
    split_log = os.getenv("SPLIT", "").lower() == "true"
    
    if split_log:
     keys = list(serializable_dict.keys())
     num_keys = len(keys)
     num_splits = math.ceil(num_keys / 3)
    
     for i in range(num_splits):
         split_keys = keys[i * 3: (i + 1) * 3]
         split_dict = {key: serializable_dict[key] for key in split_keys}

         split_filename = os.path.join(output_dir, f"updated_25_4mC_split_{i + 1}.json")
         with open(split_filename, 'w') as f:
             json.dump(split_dict, f)
         print(f"Saved {split_filename} with {len(split_dict)} DataFrames.")
    
    else:
        json_filename = os.path.join(output_dir, 'updated_25_4mC.json')
        with open(json_filename, 'w') as f:
            json.dump(serializable_dict, f)
        print("Saved updated_25_4mC.json without splitting.")

if __name__ == "__main__":
    main()
