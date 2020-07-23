"""
This plot was developed in order to plot the length (in base pairs) of crAssphage contigs before and after the filtering, using the pyplot package
"""

from Bio import SeqIO

from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = (10, 5)

contigs = SeqIO.index('all_crassphage_clean.fna', 'fasta')

contiglens = [len(contigs[contig].seq) for contig in contigs]

with open('good_contigs_ids', 'r') as f:
    good_contigs = [line.strip() for line in f.readlines()]

good_contig_lens = [len(contigs[contig].seq) for contig in contigs if contigs[contig].id in good_contigs]

plt.violinplot([good_contig_lens, contiglens], vert=False, showmedians=True)
plt.yticks((1, 2), ("viable contigs", "all contigs"))
plt.xlabel("contig length (bases)")
