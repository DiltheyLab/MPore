#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 19:16:48 2024

@author: azlannisar
"""

import pandas as pd
from Bio import SeqIO
#from intervaltree import Interval, IntervalTree


#Load Mod_Pos Data
Data_path = '/Users/azlannisar/Snakemake_Test/Sample_DF_detailed_SP10291_5mC.csv'
Mod_Data = pd.read_csv(Data_path, sep=',')

Mod_Data_plus= Mod_Data[Mod_Data['Strand']== '+' ]

fasta_file = "/Users/azlannisar/fasta_data/complete_genome/SP10291.fasta"
tsv_file = "/Users/azlannisar/TSV_REBASE_data.tsv"

with open(fasta_file, "r") as file:
    record = SeqIO.read(file, "fasta")
    sequence = str(record.seq)

#Find C Positions
c_positions = [i for i, base in enumerate(sequence) if base == 'C']
c_positions_df = pd.DataFrame(c_positions, columns=['C_Position'])
filtered_c_positions = c_positions_df[c_positions_df['C_Position'].isin(Mod_Data_plus['Mod_Pos'])]

#Initialize DF with Positions
df = pd.DataFrame(filtered_c_positions, columns=['C_Position'])

#Read tsv with Mod_Motif and Enzymes 
tsv_data = pd.read_csv(tsv_file, sep='\t')
mask = tsv_data['Motif'].str.contains('\,', case=False, na=False)
positions_with_hyphen_question = tsv_data[mask]

results_list = []

#Iterate over each Position
for c_position in df['C_Position']:
    window_start = max(0, c_position - 10)
    window_end = min(len(sequence), c_position + 10)
    
    
    for index, row in tsv_data.iterrows():
        motif = row['Mod_Motif']
        enzyme = row['Enzyme']
        
        
        motif_length = len(motif)
        for pos in range(window_start, window_end - motif_length + 1):
            if sequence[pos:pos + motif_length] == motif and c_position in range(pos, pos + motif_length):
                results_list.append({
                    'C_Position': c_position,
                    'Mod_Motif': motif,
                    'Enzyme': enzyme
                })

results_df = pd.DataFrame(results_list)

