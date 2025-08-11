#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 19:12:41 2024

@author: azlannisar
"""

import argparse
import os
import pandas as pd
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
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
        for pos, contig_number in positions:  # Unpack position and contig number
            mod_pos = pos  # Initialize mod_pos to the start position
            marked = False  # Flag to track if 'C' has been marked
            motif_with_markings = motif
            for i, char in enumerate(motif):
                if char == 'C' and not marked:
                    motif_with_markings = motif[:i] + '"' + char + '"' + motif[i+1:]
                    mod_pos += i  # Adjust mod_pos to reflect the genomic position of the marked 'C'
                    marked = True
                    break
            motif_dfs.append({
                'Motif_Start': pos,
                'Motif': motif_with_markings,
                'Mod_Pos': mod_pos,
                'Strand': strand,
                'Contig': contig_number  # Save only the numerical part
            })

            # Handle multiple 'C' in the motif
            if motif.count('C') > 1:
                next_c_index = motif.find('C', i+1)  # Find next 'C' after the current index i
                while next_c_index != -1:
                    motif_with_markings = motif[:next_c_index] + '"C"' + motif[next_c_index+1:]
                    mod_pos = pos + next_c_index  # Adjust mod_pos for the marked 'C'
                    motif_dfs.append({
                        'Motif_Start': pos,
                        'Motif': motif_with_markings,
                        'Mod_Pos': mod_pos,
                        'Strand': strand,
                        'Contig': contig_number  # Save only the numerical part
                    })
                    next_c_index = motif.find('C', next_c_index+1)  # Find next 'C'
        dfs_list.append(pd.DataFrame(motif_dfs))

    if strand == '-':
        for df in dfs_list:
            df['Motif'] = df['Motif'].apply(lambda x: x[::-1])

    return dfs_list


def save_large_dataframe(df, base_file_path, excel=False, max_entries=1_000_000):
    """
    Save a large DataFrame to multiple files if it exceeds max_entries.

    Args:
        df (pd.DataFrame): DataFrame to be saved.
        base_file_path (str): Base file path (without or with extension) to save files.
        excel (bool): Whether to save as Excel files. If False, saves as CSV files.
        max_entries (int): Maximum number of entries per file.
    """
    printed_files = set()  # Track printed files
    printed_ranges = set()  # Track printed row ranges for each file

    # Ensure the base directory exists
    base_dir = os.path.dirname(base_file_path)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    # Determine correct file extension and check if already present
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

            # Print once for each unique file and range saved
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


#Parse arguments to make thre Script run in command prompt with Input file_path being path/to/Fasta_ref
#bed_input is Path/to/Bed File which includes the results from the MM and ML Tags from the Bam File
#and output_dir being the the dir in which the results should be saved 
def parse_arguments():
    parser = argparse.ArgumentParser(description='MikrobAnalysis_5mC script')
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
    #csv_file_path = "/Users/azlannisar/Desktop/mount/Data_Test2_updated.csv"
    data_df = pd.read_csv(csv_file_path, names=["File_name", "Reference_path", "pod5_path","Bam_data","Bed_data"], skiprows=1)
    #data_df = data_df[data_df["File_name"].isin(["barcode01", "barcode02"])]
    #data_df['Reference_path'] = data_df['Reference_path'].str.replace(
    #"/home/azlan/VRE_Data/ref/Reference/",
    #"/Users/azlannisar/Desktop/mount2/ref/Reference/")
    #data_df['Bed_data'] = data_df['Bed_data'].str.replace(
    #"/home/azlan/VRE_Data/Output/",
    #"/Users/azlannisar/Desktop/mount/")
    
    #Adding the motifs from rebase which are only in the homology search and defined by user 
    
    include_rebase_motifs = os.getenv("INCLUDE_REBASE_MOTIFS", "False").lower() == "true"
    motifs_file_path = args.Motif_list
    #motifs_file_path = "/Users/azlannisar/Desktop/mount2/Motifs_6mA.txt"
    motifs = []
    if motifs_file_path:
        with open(motifs_file_path, "r") as file:
            motifs = [line.strip() for line in file if line.strip()]
    if include_rebase_motifs:
        #Include Motifs from REBASE If needed not all 
        mtase_file = args.Mtase_File
        #mtase_file = "/Users/azlannisar/Desktop/mount/Mtase_presence_e_25_values.csv"
        mtase_df = pd.read_csv(mtase_file)
        enzyme_names = mtase_df.columns[1:].tolist() 
        
        tsv_file = args.REBASE_Motifs
        #tsv_file = "/Users/azlannisar/Desktop/mount2/TSV_Enzyme.csv"
        tsv_df = pd.read_csv(tsv_file, sep="\t")
        tsv_df = tsv_df[tsv_df['MethylationType'] == "5mC"]
        
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
        
        filtered_motifs = [
        motif for motif in filtered_motifs 
        if re.match(r'^[A-Za-z]+$', motif) and len(motif) > 3]
        
        motifs = list(set(motifs + filtered_motifs))
        
    reversed_motifs = [motif[::-1] for motif in motifs]
    
    #Check Block for csv file (directories and reference files..)
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
        filename = row['File_name']  # Get the File_name
        file_path = row['Reference_path']  # Get the Reference_path
        fasta_dict = {}  # Create a new dictionary for each row (File_name)
        
        # Initialize dictionaries to accumulate motif matches across all sequences for each motif
        matches_forward_all = {}
        matches_reverse_all = {}
        
        with open(file_path, "r") as fasta_file:
            for fasta in SeqIO.parse(fasta_file, "fasta"):
                fasta_dict[fasta.id] = str(fasta.seq)
        
                # Get contig sequence and reverse complement
                contig_sequence = fasta_dict[fasta.id].upper()
                contig_sequence_rev = str(Seq(contig_sequence).complement()).upper()
        
                # Extract numerical part of the contig ID
                contig_number = fasta.id  # Extract number
        
                # Search for motifs
                matches_forward = search_motifs(contig_sequence, motifs)
                matches_reverse = search_motifs(contig_sequence_rev, reversed_motifs)
        
                # Save matches with contig number
                for motif, positions in matches_forward.items():
                    for pos in positions:
                        if motif not in matches_forward_all:
                            matches_forward_all[motif] = []
                        matches_forward_all[motif].append((pos, contig_number))  # Save position and contig number
        
                for motif, positions in matches_reverse.items():
                    for pos in positions:
                        if motif not in matches_reverse_all:
                            matches_reverse_all[motif] = []
                        matches_reverse_all[motif].append((pos, contig_number)) 
    
        # After processing all sequences, convert sets back to lists
        matches_forward_all = {motif: sorted(list(positions)) for motif, positions in matches_forward_all.items()}
        matches_reverse_all = {motif: sorted(list(positions)) for motif, positions in matches_reverse_all.items()}
       # matches_forward_all = {k: v for k, v in matches_forward_all.items() if len(v) >= 20 and len(k) >= 3}
       # matches_reverse_all = {k: v for k, v in matches_reverse_all.items() if len(v) >= 20 and len(k) >= 3}
        
        # Generate combined dataframes for forward and reverse matches
        dfs = generate_motif_dataframes(matches_forward_all, '+')
        dfs_rev = generate_motif_dataframes(matches_reverse_all, '-')
    
    
    #for index, row in data_df.iterrows():
        #filename = row['File_name']  # Get the File_name
        #file_path = row['Reference_path']  # Get the Reference_path
        #fasta_dict = {}
        #matches_forward_all = {}
        #matches_reverse_all = {}
        # Open and read the FASTA file
        #with open(file_path, "r") as fasta_file:
            #for fasta in SeqIO.parse(fasta_file, "fasta"):
             #   fasta_dict[fasta.id] = str(fasta.seq)
        #contig_33_sequence = fasta_dict.get(fasta.id, '')
        #contig_33_sequence_rev = str(Seq(contig_33_sequence).complement())
        #if contig_33_sequence.islower():
         #   contig_33_sequence = contig_33_sequence.upper()
        #if contig_33_sequence_rev.islower():
         #   contig_33_sequence_rev = contig_33_sequence_rev.upper()
    
        
        #matches_forward = search_motifs(contig_33_sequence, motifs)
        #matches_reverse = search_motifs(contig_33_sequence_rev, reversed_motifs)
        #matches_forward = {k: v for k, v in matches_forward.items() if len(v) >= 20 and len(k) >= 3}
        #matches_reverse = {k: v for k, v in matches_reverse.items() if len(v) >= 20 and len(k) >= 3}
        #dfs = generate_motif_dataframes(matches_forward,'+')
        #dfs_rev = generate_motif_dataframes(matches_reverse, '-')
        
        #for df in dfs:
            #df["Mod_Pos"] = df["Mod_Pos"] +1
        #for df in dfs_rev:
            #df["Mod_Pos"] = df ["Mod_Pos"] +1
    
        #Merge Dfs from forward and reverse matches according to their Motif. e.g CCAGG on forward 
        #has the match GGACC on reverse. 
        merged_dfs = []

        for i in range(len(dfs)):
            merged_df = pd.concat([dfs[i], dfs_rev[i]], axis=0) 
            merged_dfs.append(merged_df)
    
    
        #motif_ccwgg = ['CCAGG', 'CCTGG']
        #if any(motif in args.Motif_list for motif in motif_ccwgg):
            #ccwgg_df = filter_and_concat(merged_dfs, motif_ccwgg)
            #ccwgg_df['motif'] = 'CCWGG'
            #ccwgg_df.loc[ccwgg_df['strand'] == '-', 'motif'] = ccwgg_df.loc[ccwgg_df['strand'] == '-', 'motif'].apply(lambda x: x[::-1])
            #merged_dfs.append(ccwgg_df)
        
        #motif_gcngc = ['GCAGC', 'GCCGC', 'GCGGC', 'GCTGC']
        #if any(motif in args.Motif_list for motif in motif_gcngc):
            #gcngc_df = filter_and_concat(merged_dfs, motif_gcngc)
            #gcngc_df['motif'] = 'GCNGC'
            #gcngc_df.loc[gcngc_df['strand'] == '-', 'motif'] = gcngc_df.loc[gcngc_df['strand'] == '-', 'motif'].apply(lambda x: x[::-1])
            #gcngc_df = gcngc_df.drop(columns=['Mod_Pos3'])
            #merged_dfs.append(gcngc_df)
        
        #motif_ggncc = ['GGACC', 'GGCCC', 'GGGCC', 'GGTCC']
        #if any(motif in args.Motif_list for motif in motif_ggncc):
            #ggncc_df = filter_and_concat(merged_dfs, motif_ggncc)
            #ggncc_df['motif'] = 'GGNCC'
            #ggncc_df.loc[ggncc_df['strand'] == '-', 'motif'] = ggncc_df.loc[ggncc_df['strand'] == '-', 'motif'].apply(lambda x: x[::-1])
            #ggncc_df = ggncc_df.drop(columns=['Mod_Pos3'])
            #merged_dfs.append(ggncc_df)
        
        #motif_gantc = ['GAATC', 'GACTC', 'GAGTC', 'GATTC']
        #if any(motif in args.Motif_list for motif in motif_gantc):
            #gantc_df = filter_and_concat(merged_dfs, motif_gantc)
            #gantc_df['motif'] = 'GANTC'
            #gantc_df.loc[gantc_df['strand'] == '-', 'motif'] = gantc_df.loc[gantc_df['strand'] == '-', 'motif'].apply(lambda x: x[::-1])
            #gantc_df = gantc_df.drop(columns=['Mod_Pos2'])
            #merged_dfs.append(gantc_df)    
    
        #Read Bed File containing the Cols 'start', 'end' and 'score' at least. Since they are
        #used downstream in the analysis
        # Read the TSV file into a DataFrame
        
        chunk_size = 100000  # Adjust as needed
        Barcode06_Wig_5mC = pd.DataFrame()
        
        # Process in chunks
        for chunk in pd.read_csv(row['Bed_data'], sep='\t', header=None, chunksize=chunk_size):
            chunk = chunk.reset_index(drop=True)
            
            # Check if column 9 contains strings with spaces and split if needed
            if chunk[9].apply(lambda x: isinstance(x, str)).any():
                if chunk[9].str.contains(' ').any():
                    split_values = chunk[9].str.split(expand=True)
                    chunk.drop(columns=[9], inplace=True)
                    split_values.columns = [f"split_{i}" for i in range(split_values.shape[1])]
                    chunk = pd.concat([chunk, split_values], axis=1)
            
            # Check if the dataframe has string headers; if not, add them
            headers = ['sample', 'start', 'end', 'name', 'score', 'strand', 'tstart', 'tend', 'color', 
                       'coverage', 'frequence', 'modified', 'canonical', 'other_mod', 'delted', 
                       'failed', 'substitution', 'No_call']
            
            if not chunk.iloc[0].apply(lambda x: isinstance(x, str)).all():
                chunk.columns = headers[:len(chunk.columns)]
            
            # Filter rows where 'name' is 'm'
            chunk = chunk[chunk['name'] == 'm']
            
            # Split into positive and negative strands
            positive_strand_data = chunk[chunk['strand'] == '+']
            negative_strand_data = chunk[chunk['strand'] == '-']
            
            # Create dataframes for downstream analysis
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
            
            # Concatenate results
            Barcode06_Wig_5mC = pd.concat([Barcode06_Wig_5mC, positive_df, negative_df], ignore_index=True)
    
        #Add methylation scores to the merged dataframes occording to the positions
        #Important to see if it is 0-indexed or 1-indexed else adjust the code with
        #Space in Ref Seq string or adjust the Positions with df.loc[:,column_to_adjust] -1 
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
    
        #Filter Positions in which no methylation score was detected. Should be almost 0 else check 
        #Base quality and or reference sequence
        filtered_dfs = []
        for df in merged_dfs:
            filtered_df = df[df['Score1'] != '/']
            filtered_dfs.append(filtered_df)
    
        #Add dfs from Filtered_list to Boxplot_list. Also the Motifs are being seperated and used as
        #Key value. E.g. the Dataframe with the Motifs CCAGG and its matching GGACC or CCAGG- on reverse
        #Are now saved seperately and with their Key Value is set occording to the Motifs.
        #Also a pointer is being added to the Key values as integer which shows which dfs are matches
        #occurding to forward and reverse
        Boxplot_dfs = {}

        for df in filtered_dfs:
            for motif in df['Motif'].unique():
                for strand in df['Strand'].unique():
                    for contig in df['Contig'].unique():
                        subset_df = df[(df['Motif'] == motif) & (df['Strand'] == strand)]
                        key = f"{motif}({strand})"
                        if key not in Boxplot_dfs:
                            Boxplot_dfs[key] = subset_df
                        else:
                            Boxplot_dfs[key] = pd.concat([Boxplot_dfs[key], subset_df])
        for key, df in Boxplot_dfs.items():
            Boxplot_dfs[key] = df.drop_duplicates(subset=['Mod_Pos', 'Motif', 'Strand', 'Contig'])                    
    
        #Get the filename from parse file_path and remove the extension e.g. SP10291_RIIPC_n.fasta
        #becomes SP10291_RIIPC_n it is advised to see name the reference files occordingly so the 
        #output names are distinct 
        list_name = filename
        dir_path = os.getenv("OUTPUT_DIR")
        #dir_path= "/Users/azlannisar/"
        #Create Output_Dir if it does not exist
        os.makedirs(dir_path, exist_ok=True)
    
        #Export the List to output_path
        list_output_path = os.path.join(dir_path, f"{list_name}_5mC.pkl")
        with open(list_output_path, 'wb') as f:
            pickle.dump(Boxplot_dfs, f)
        
    
        #Initialize the Sample_DF with Motif and matching Key. The Motifs are from Boxplot_dfs
        #each Motif is being added matching to their Key. If a Motif has two or more Cs the keys
        #are added more often e.g. CC"A"GG and CC"A"GG are both added under the key CCAGG.
        #The number of Sites for each Motif in the Ref Sequence is also added to the Dataframe. 
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
        rows_to_keep = concatenated_df.apply(lambda row: len(row) <= 6, axis=1)
        concatenated_df = concatenated_df[rows_to_keep]
        concatenated_df.reset_index(drop=True, inplace=True)
                
        list_output_path_2 = os.path.join(dir_path, f"{list_name}_5mC_merged.pkl")
        
        with open(list_output_path_2, 'wb') as f:
            pickle.dump(merged_dfs_copy, f)
        
        Sample_DF = pd.DataFrame(Sample_DF) 
        concatenated_df = pd.DataFrame(concatenated_df)  
        
        #Generate the Name of the Files being saved and save them in their output_dir. 
        xlsx_file_path = os.path.join(dir_path, f"Sample_DF_{list_name}_5mC.xlsx")
        csv_file_path = os.path.join(dir_path, f"Sample_DF_{list_name}_5mC.csv")
        
        xlsx_file_path2 = os.path.join(dir_path, f"Sample_DF_detailed_{list_name}_5mC")
        csv_file_path2 = os.path.join(dir_path, f"Sample_DF_detailed_{list_name}_5mC")
        
        csv_file_path3 = os.path.join(dir_path, f"Sample_DF_detailed_{list_name}_5mC.csv")
        
        
        save_large_dataframe(Sample_DF, xlsx_file_path, excel=True)
        save_large_dataframe(Sample_DF, csv_file_path, excel=False)
        save_large_dataframe(concatenated_df, xlsx_file_path2, excel=True)
        save_large_dataframe(concatenated_df, csv_file_path2, excel=False)
        
        concatenated_df.to_csv(csv_file_path3, index=False)
    
if __name__ == "__main__":
    main()
