"""
A simple script for removing duplicate sequences with different headers from a FASTA file.
It makes use of the Biopython toolbox.
"""

from Bio import SeqIO
from itertools import combinations

all_cr_raw = list(SeqIO.parse('all_crassphage.fna', 'fasta'))

duplicate_positions = []

for i, j in combinations(range(len(all_cr_raw)), 2):
    if all_cr_raw[i].seq == all_cr_raw[j].seq:
        duplicate_positions.append(j)

all_cr_clean = [all_cr_raw[i] for i in range(len(all_cr_raw)) if i not in duplicate_positions]
        
SeqIO.write(all_cr_clean, 'all_crassphage_clean.fna', 'fasta')
