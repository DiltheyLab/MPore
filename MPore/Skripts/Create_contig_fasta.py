#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 16:55:36 2024

@author: azlannisar
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

input_fasta_path = "/Users/azlannisar/fasta_data/complete_genome/SP10291.fasta"
output_fasta_path = "/Users/azlannisar/fasta_data/complete_genome/SP10291_two_contigs.fasta"

records = list(SeqIO.parse(input_fasta_path, "fasta"))

if len(records) != 1:
    raise ValueError("The original FASTA file should contain exactly one sequence.")

original_seq = str(records[0].seq)

half_length = len(original_seq) // 2
first_half = original_seq[:half_length]
second_half = original_seq[half_length:]
contig1 = SeqRecord(Seq(first_half), id="contig1", description="First half of original sequence")
contig2 = SeqRecord(Seq(second_half), id="contig2", description="Second half of original sequence")
with open(output_fasta_path, "w") as output_handle:
    SeqIO.write([contig1, contig2], output_handle, "fasta")