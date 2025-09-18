#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 19:15:17 2024

@author: azlannisar
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
import numpy as np


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

def iupac_to_regex_motif(motif):
    regex_motif = ''
    for base in motif:
        regex_motif += iupac_to_regex.get(base, base)
    return regex_motif

#Function to search Motifs in the Fasta_Ref Sequence String
def search_motifs(sequence, motifs):
    matches = {}
    for motif in motifs:
        regex_motif = iupac_to_regex_motif(motif)
        motif_positions = [m.start() for m in re.finditer(regex_motif, sequence)]
        matches[motif] = motif_positions
    return matches

#Function to generate Dataframe from the matches found in Fasta_Ref Sequence String,
#It also applies Mod_Pos to see which Position in the reference Genome is modified 
def generate_motif_dataframes(matches, strand):
    dfs_list = []
    for motif, positions in matches.items():
        motif_dfs = []
        for pos, contig_number in positions:  
            mod_pos = pos  
            marked = False
            motif_with_markings = motif
            for i, char in enumerate(motif):
                if char == 'A' and not marked:
                    motif_with_markings = motif[:i] + '"' + char + '"' + motif[i+1:]
                    mod_pos += i  
                    marked = True
                    break
            motif_dfs.append({
                'Motif_Start': pos,
                'Motif': motif_with_markings,
                'Mod_Pos': mod_pos,
                'Strand': strand,
                'Contig': contig_number  
            })

            if motif.count('A') > 1:
                next_c_index = motif.find('A', i+1)  
                while next_c_index != -1:
                    motif_with_markings = motif[:next_c_index] + '"A"' + motif[next_c_index+1:]
                    mod_pos = pos + next_c_index  
                    motif_dfs.append({
                        'Motif_Start': pos,
                        'Motif': motif_with_markings,
                        'Mod_Pos': mod_pos,
                        'Strand': strand,
                        'Contig': contig_number  
                    })
                    next_c_index = motif.find('A', next_c_index+1)  
        dfs_list.append(pd.DataFrame(motif_dfs))

    if strand == '-':
        for df in dfs_list:
            df['Motif'] = df['Motif'].apply(lambda x: x[::-1])

    return dfs_list

def save_large_dataframe(df, base_file_path, excel=False, max_entries=1_000_000):
    printed_files = set()  
    printed_ranges = set()  

    base_dir = os.path.dirname(base_file_path)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    file_extension = ".xlsx" if excel else ".csv"
    if not base_file_path.endswith(file_extension):
        base_file_path = base_file_path + file_extension

    total_rows = len(df)
    if total_rows > max_entries:
        num_files = int(np.ceil(total_rows / max_entries))
        for i in range(num_files):
            start_row = i * max_entries
            end_row = min((i + 1) * max_entries, total_rows)
            sub_df = df.iloc[start_row:end_row]
            sub_file_path = f"{base_file_path.rsplit('.', 1)[0]}_{i + 1}{file_extension}"

            if excel:
                sub_df.to_excel(sub_file_path, index=False)
            else:
                sub_df.to_csv(sub_file_path, index=False)


            if sub_file_path not in printed_files or (start_row, end_row) not in printed_ranges:
                print(f"Saved {sub_file_path} with rows from {start_row} to {end_row - 1}")
                printed_files.add(sub_file_path)
                printed_ranges.add((start_row, end_row))
                
    else:
        if excel:
            df.to_excel(base_file_path, index=False)
        else:
            df.to_csv(base_file_path, index=False)
        if base_file_path not in printed_files:
            print(f"Saved {base_file_path} with all rows")
            printed_files.add(base_file_path)
        
   
        

def parse_arguments():
    parser = argparse.ArgumentParser(description='MikrobAnalysis script')
    parser.add_argument('csv_list', type=str, help='Path to .csv file')
    parser.add_argument('output_dir', type=str, help='Output directory path')
    parser.add_argument('Motif_list', type=str, help='List_of_Motifs')
    parser.add_argument('Mtase_File', type=str, help="MTase_File_from_Blast")
    parser.add_argument('TSV_File', type=str, help="TSV_File_from_REBASE")
    parser.add_argument('REBASE_Motifs', type=str, help="REBASE_Motifs_TSV")
    return parser.parse_args()

#Main Function 
def main():
    
    args = parse_arguments()
    csv_file_path= args.csv_list
    data_df = pd.read_csv(csv_file_path, names=["File_name", "Reference_path", "pod5_path","Bam_data","Bed_data"], skiprows=1)    
    include_rebase_motifs = os.getenv("INCLUDE_REBASE_MOTIFS", "true").lower() == "true"
    print(f"Include REBASE Motifs: {include_rebase_motifs}")
    motifs_file_path = args.Motif_list
    print(f"Motifs file path: {motifs_file_path}")
    motifs = []
    if motifs_file_path:
        with open(motifs_file_path, "r") as file:
            motifs = [line.strip() for line in file if line.strip()]
    print(f"Loaded motifs: {motifs}")
    if include_rebase_motifs:
        mtase_file = args.Mtase_File
        mtase_df = pd.read_csv(mtase_file)
        enzyme_names = mtase_df.columns[1:].tolist()
        print(f"Enzyme Names: {enzyme_names}") 
        
        tsv_file = args.REBASE_Motifs
        tsv_df = pd.read_csv(tsv_file, sep="\t")
        tsv_df = tsv_df[tsv_df['MethylationType'] == "6mA"]
        print(tsv_df.head())
        filtered_motifs = []
        
        for enzyme in enzyme_names:
            enzyme_df = tsv_df[tsv_df['Enzyme'] == enzyme]
            for i, row in enzyme_df.iterrows():
                if '?' in str(row['Position']) or '-' in str(row['Strand']):
                    motif = row['Motif']
                    complement_motif = str(Seq(motif).complement())[::-1]
                    filtered_motifs.append(motif)
                    filtered_motifs.append(complement_motif)
                else:
                    filtered_motifs.append(row['Motif'])
        print(f"Filtered Motifs: {filtered_motifs}")
        filtered_motifs = [
        motif for motif in filtered_motifs 
        if re.match(r'^[A-Za-z]+$', motif) and len(motif) > 3]
        print(f"Motifs before adding REBASE: {len(motifs)}")
        motifs = list(set(motifs + filtered_motifs))
        print(f"Motifs after adding REBASE: {len(motifs)}")
    reversed_motifs = [motif[::-1] for motif in motifs]
    

    for index, row in data_df.iterrows():
        reference_path= row['Reference_path']
        bed_data_path= row['Bed_data']
        if not os.path.exists(reference_path):
            print(f"Error: Reference directory not found. Index: {index}, Reference Path: {reference_path}")
        if not os.path.exists(bed_data_path):
            print(f"Error: Bed data directory not found. Index: {index}, Bed Data Path: {bed_data_path}")
        if not reference_path.lower().endswith(('.fasta','.fa')):
            print(f"Error: Reference file should be a .fasta or .fa file. Index: {index}")
    
    
    for index, row in data_df.iterrows():
        filename = row['File_name']  
        file_path = row['Reference_path']  
        fasta_dict = {}  
        
        matches_forward_all = {}
        matches_reverse_all = {}
        
        with open(file_path, "r") as fasta_file:
            for fasta in SeqIO.parse(fasta_file, "fasta"):
                contig_id = fasta.id

                contig_sequence = str(fasta.seq).upper()
                contig_sequence_rev = str(Seq(contig_sequence).complement()).upper()
                contig_number = fasta.id  
                matches_forward = search_motifs(contig_sequence, motifs)
                matches_reverse = search_motifs(contig_sequence_rev, reversed_motifs)
                for motif, positions in matches_forward.items():
                    for pos in positions:
                        if motif not in matches_forward_all:
                            matches_forward_all[motif] = []
                        matches_forward_all[motif].append((pos, contig_number)) 
        
                for motif, positions in matches_reverse.items():
                    for pos in positions:
                        if motif not in matches_reverse_all:
                            matches_reverse_all[motif] = []
                        matches_reverse_all[motif].append((pos, contig_number)) 
                        
        

        matches_forward_all = {motif: sorted(list(positions)) for motif, positions in matches_forward_all.items()}
        matches_reverse_all = {motif: sorted(list(positions)) for motif, positions in matches_reverse_all.items()}
        dfs = generate_motif_dataframes(matches_forward_all, '+')
        dfs_rev = generate_motif_dataframes(matches_reverse_all, '-')
        merged_dfs = []

        for i in range(len(dfs)):
            merged_df = pd.concat([dfs[i], dfs_rev[i]], axis=0)
            merged_dfs.append(merged_df)

        chunk_size = 100000  
        Barcode06_Wig_5mC = pd.DataFrame()
        

        for chunk in pd.read_csv(row['Bed_data'], sep='\t', header=None, chunksize=chunk_size):
            chunk = chunk.reset_index(drop=True)
            
            if chunk[9].apply(lambda x: isinstance(x, str)).any():
                if chunk[9].str.contains(' ').any():
                    split_values = chunk[9].str.split(expand=True)
                    chunk.drop(columns=[9], inplace=True)
                    split_values.columns = [f"split_{i}" for i in range(split_values.shape[1])]
                    chunk = pd.concat([chunk, split_values], axis=1)
            

            headers = ['sample', 'start', 'end', 'name', 'score', 'strand', 'tstart', 'tend', 'color', 
                       'coverage', 'frequence', 'modified', 'canonical', 'other_mod', 'delted', 
                       'failed', 'substitution', 'No_call']
            
            if not chunk.iloc[0].apply(lambda x: isinstance(x, str)).all():
                chunk.columns = headers[:len(chunk.columns)]

            chunk = chunk[chunk['name'] == 'a']
            
            positive_strand_data = chunk[chunk['strand'] == '+']
            negative_strand_data = chunk[chunk['strand'] == '-']
            
            positive_df = pd.DataFrame({
                'Mod_Pos': positive_strand_data['start'],
                'Meth_Score_5mC': positive_strand_data['frequence'],
                'Strand': '+',
                'contig': positive_strand_data['sample']
            })
            negative_df = pd.DataFrame({
                'Mod_Pos': negative_strand_data['start'],
                'Meth_Score_5mC': negative_strand_data['frequence'],
                'Strand': '-',
                'contig': negative_strand_data['sample']
            })
            Barcode06_Wig_5mC = pd.concat([Barcode06_Wig_5mC, positive_df, negative_df], ignore_index=True)

        score_lookup = {}
        for _, row in Barcode06_Wig_5mC.iterrows():
            key = (row['Mod_Pos'], row['Strand'], row['contig'])
            score_lookup[key] = row['Meth_Score_5mC']
        for df in merged_dfs:
            scores1 = []
    
            for _, row in df.iterrows():
                key1 = (row['Mod_Pos'], row['Strand'], row['Contig'])
                score1 = score_lookup.get(key1, '/')
                scores1.append(score1)
        
            df['Score1'] = scores1


        filtered_dfs = []
        for df in merged_dfs:
            filtered_df = df[df['Score1'] != '/']
            filtered_dfs.append(filtered_df)
    
        Boxplot_dfs = {}

        for df in filtered_dfs:
            for motif in df['Motif'].unique():
                for strand in df['Strand'].unique():
                    subset_df = df[(df['Motif'] == motif) & (df['Strand'] == strand)]
                    key = f"{motif}({strand})"
                    if key not in Boxplot_dfs:
                        Boxplot_dfs[key] = subset_df
                    else:
                        Boxplot_dfs[key] = pd.concat([Boxplot_dfs[key], subset_df])
        
        for key, df in Boxplot_dfs.items():
            Boxplot_dfs[key] = df.drop_duplicates(subset=['Mod_Pos', 'Motif', 'Strand', 'Contig'])
        list_name = filename
        dir_path = args.output_dir

        os.makedirs(dir_path, exist_ok=True)
    
        num_rows_dict = {key: len(value) for key, value in Boxplot_dfs.items()}
        
        keys = list(Boxplot_dfs.keys())
        Sample_DF = pd.DataFrame({'Key': keys})
        Sample_DF[f"Number_of_Sites_{list_name}"] = Sample_DF['Key'].map(num_rows_dict)
        for key, value in Boxplot_dfs.items():
            value.loc[:, 'Score1'] = value['Score1'].astype(float)
        mean_score_dict = {key: value['Score1'].mean() for key, value in Boxplot_dfs.items()}
        Sample_DF['Avg_Methylation'] = Sample_DF['Key'].map(mean_score_dict)

        merged_dfs_copy = merged_dfs.copy()
        concatenated_df = pd.concat(merged_dfs_copy, ignore_index=True)

        Sample_DF = pd.DataFrame(Sample_DF) 
        concatenated_df = pd.DataFrame(concatenated_df) 
        rows_to_keep = concatenated_df.apply(lambda row: len(row) <= 6, axis=1)
        concatenated_df = concatenated_df[rows_to_keep]

        xlsx_file_path = os.path.join(dir_path, f"Sample_DF_{list_name}_6mA.xlsx")
        csv_file_path = os.path.join(dir_path, f"Sample_DF_{list_name}_6mA.csv")
        
        xlsx_file_path2 = os.path.join(dir_path, f"Sample_DF_detailed_{list_name}_6mA")
        csv_file_path2 = os.path.join(dir_path, f"Sample_DF_detailed_{list_name}_6mA")
        
        csv_file_path3 = os.path.join(dir_path, f"Sample_DF_detailed_{list_name}_6mA.csv")

        
        save_large_dataframe(Sample_DF, xlsx_file_path, excel=True)
        save_large_dataframe(Sample_DF, csv_file_path, excel=False)
        save_large_dataframe(concatenated_df, xlsx_file_path2, excel=True)
        save_large_dataframe(concatenated_df, csv_file_path2, excel=False)
        
        concatenated_df.to_csv(csv_file_path3, index=False)
    
if __name__ == "__main__":
    main()
