"""
This script reads every FASTA alignment file and calculates, for each position, the information quantity.
The data is then written in a file, one for each alignment.
It makes use of the Biopython AlignIO  module
"""

from os import listdir
from math import log2
from Bio import AlignIO

amino_acids = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

for filename in [file for file in listdir() if '.msa' in file]:
    
    alignment = AlignIO.read(filename, 'fasta')
    
    with open(filename.replace('msa', 'csv'), 'w') as output:
        output.write('SUM' + '\t')
        output.write('coverage' + '\t')
        output.write('\t'.join(amino_acids))
        output.write('\n')
        for i in range(alignment.get_alignment_length()):
            output_line = []
            column = alignment[:, i]
            howmany_seqs = sum([1 for x in column if x != '-'])
            for letter in amino_acids:
                if column.count(letter) == 0:
                    output_line.append(0)
                else:    
                    rate = column.count(letter)/len(column)
                    bit = rate*log2(rate*20)
                    output_line.append(bit)
            output.write(str(sum(output_line)) + '\t')
            output.write(str(howmany_seqs) + '\t')
            output.write('\t'.join([str(x) for x in output_line]))
            output.write('\n')
