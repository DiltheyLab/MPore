#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 17:32:05 2024

@author: azlannisar
"""

import pandas as pd

# Beispiel-DataFrame
data = {
    'Enzyme': ['M.PSPGI', 'M.PSPGI2', 'M.PARA', 'M.DIFF', 'M.EASY', 'M.EXAMP'],
    'Motif': ['CCWGG', 'CCWGG', 'CCWTG', 'CCWTG', 'CCWTG', 'CCWCG'],
    'Position': ['2(5)', '?(5)', '?(5),-4(6)', '?(5),?(5)', '1(5),-4(6)', '4(5),?(5)']
}

df = pd.DataFrame(data)

# Funktion zur Extraktion der Informationen
def extract_info(position):
    positions = []
    methyl_types = []
    strands = []
    
    # Split positions by comma
    entries = position.split(',')
    
    for entry in entries:
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

# Neue Spalten initialisieren
df['NewPosition'] = None
df['MethylType'] = None
df['Strand'] = None

# Informationen extrahieren und in neue Spalten speichern
for i, row in df.iterrows():
    positions, methyl_types, strands = extract_info(row['Position'])
    df.at[i, 'NewPosition'] = ','.join(positions)
    df.at[i, 'MethylType'] = ','.join(methyl_types)
    df.at[i, 'Strand'] = ','.join(strands)

# DataFrame anzeigen
print(df)
