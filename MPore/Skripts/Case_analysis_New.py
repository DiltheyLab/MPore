#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 16:00:26 2024

@author: azlannisar
"""

import pandas as pd

# Beispiel-DataFrame
data = {
    'Enzyme': ['M.PSPGI', 'M.PSPGI2', 'M.PARA', 'M.DIFF', 'M.EASY', 'M.EXAMP'],
    'Motif': ['CCWGG', 'CCWGG', 'CCWTG', 'CCWTG', 'CCWTG', 'CCWCG'],
    'Position': ['2(5)', '?(5)', '?(5),-4(6)', '?(5),?(5)', '1(5),-4(6)', '4(5),?(5)'],
    'NewPosition': ['2', '?', '?,4', '?,?', '1,4', '4,?'],
    'MethylType': ['5mC', '5mC', '5mC,6mA', '5mC,5mC', '5mC,6mA', '5mC,5mC'],
    'Strand': ['+', '+', '+,-', '+,+', '+,-', '+,+']
}

df = pd.DataFrame(data)
df.loc[df['NewPosition'] == '?,?', 'NewPosition'] = '?'

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

# Funktion zur Verarbeitung der Daten
def process_data(df):
    result_rows = []

    for _, row in df.iterrows():
        enzyme = row['Enzyme']
        motif = row['Motif']
        positions = row['NewPosition'].split(',')
        methyl_types = row['MethylType'].split(',')
        strands = row['Strand'].split(',')

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

        # 5 Case: Eine bekannte und eine unbekannte Position auf demselben Strand mit identischem Methyltyp
        elif len(positions) == 2 and '?' in positions and len(set(strands)) == 1    :
            unknown_index = positions.index('?')
            known_index = 1 - unknown_index
            motif_list = list(motif)
            if positions[known_index].isdigit():
                motif_list[int(positions[known_index]) - 1] = 'X'  # Ersetzen des bekannten C durch X im Motif
                new_motif = ''.join(motif_list)
                result_rows.append([enzyme, f"{enzyme}_C{positions[known_index]}", motif, methyl_types[known_index], strands[known_index], positions[known_index]])
                result_rows.extend(generate_combinations(new_motif, methyl_types[unknown_index], enzyme, strands[unknown_index]))

    # Erstellung eines neuen DataFrames mit den Ergebnissen
    result_df = pd.DataFrame(result_rows, columns=['Enzyme', 'EnzymeUniqueMethylationActivity', 'Motif', 'MethylationType', 'Strand', 'Position'])
    return result_df

# Verarbeiten der Daten und Ausgabe des neuen DataFrames
result_df = process_data(df)
print(result_df)
