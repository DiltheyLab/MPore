#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 18:25:55 2024

@author: azlannisar
"""

import pandas as pd
from Bio.Seq import Seq

def reverse_complement_motif(motif):
    return str(Seq(motif).reverse_complement())

tsv_file_path = "/Users/azlannisar/Type2_Datframe.tsv"
data = pd.read_csv(tsv_file_path, sep='\t')

entries_to_remove = [
    '(10/12)CGANNNNNNTGC(12/10)',
    '(10/15)ACNNNNGTAYC(12/7)',
    '(11/13)CAANNNNNGTGG(12/10)',
    '(11/13)GACNNNNRTGA(12/10)',
    '(12/14)GACNNNNNTGA(13/11)',
    '(7/12)GAACNNNNNCTC(13/8)',
    '(7/12)GAACNNNNNNTCC(12/7)',
    '(7/13)GAYNNNNNRTC(14/9)',
    '(8/13)CACNNNNNNTCC(12/7)',
    '(8/14)CAYNNNNNRTG(14/8)',
    '(8/14)CCANNNNNNGT(15/9)',
    '(9/11)TGRYCA(11/9)',
    ')01/21(CCACNNNNNTTG)31/11(',
    ')01/21(GCANNNNNNTCG)21/01(',
    ')7/21(GGANNNNNNGTG)31/8(',
    ')7/21(GGANNNNNNGTTC)21/7(',
    ')7/21(GRTACNNNNGT)51/01(',
    ')7/9(CTGAG',
    ')8/01(CTCCTC',
    ')8/31(GAGNNNNNGTTC)21/7('
]

data = data[~data['Motif'].isin(entries_to_remove)]

data['Position'] = data['Position'].astype(str)

# Identify rows where Position column contains ',-' and duplicate those rows
indices_to_duplicate = data[data['Position'].str.contains(',-')].index
duplicated_rows = data.loc[indices_to_duplicate].copy()

# Add a suffix to the Enzyme column of the duplicated rows to differentiate them
duplicated_rows['Enzyme'] = duplicated_rows['Enzyme'] + '_reverse'

duplicated_rows['Mod_Motif'] = duplicated_rows['Motif'].apply(reverse_complement_motif)
data['Mod_Motif'] = data['Motif']


df = pd.concat([data, duplicated_rows], ignore_index=True)
output_csv_path = "/Users/azlannisar/TSV_REBASE_data.tsv"
df.to_csv(output_csv_path, sep='\t', index=False)

# Print the resulting DataFrame (optional)
print(df)

