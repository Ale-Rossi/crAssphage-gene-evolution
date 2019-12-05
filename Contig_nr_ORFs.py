"""
This script was used for plotting the number of predicted ORFs in crAss-like contigs before and after the filtering.
It makes use of the pyplot package
"""

from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = (10, 5)

def contigname(seqname):
    return '_'.join(seqname.split('_')[:-1])
def seqnumber(seqname):
    return int(seqname.split('_')[-1])

with open('all_crassphage_clean.fna', 'r') as f:
    all_clean_contigs = [line.strip().replace('>', '').split(' ')[0] for line in f.readlines() if '>' in line]

with open('all_crassphage_clean.pred.fna', 'r') as f:
    all_orfs = [line.strip().replace('>', '').split(' ')[0] for line in f.readlines() if '>' in line]

contig_orfs = {contig: [] for contig in all_clean_contigs}
for orf in all_orfs:
    contig_orfs[contigname(orf)].append(orf)

contig_lens = [len(contig_orfs[contig]) for contig in contig_orfs]

with open('good_contigs_ids', 'r') as f:
    good_contigs = [line.strip() for line in f.readlines()]

good_contig_lens = [len(contig_orfs[contig]) for contig in good_contigs]

plt.violinplot([good_contig_lens, contig_lens], vert=False, showmedians=True)
plt.yticks((1, 2), ("viable contigs", "all contigs"))
plt.xlabel("n° of ORFs")
