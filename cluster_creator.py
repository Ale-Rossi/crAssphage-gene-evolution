"""
A script used to create the alignment units from the PSI-BLAST tabular output.
It makes use of the Biopython SeqIO module.
"""

from Bio import SeqIO

# The PSI-BLAST output file is read; lines that do not contain matches are not considered
with open('/home/alessandro/crassphage_def/all_crassphage_clean.blast') as f:
    blast_output = [line.strip().split('\t') for line in f.readlines() if len(line.strip().split('\t')) > 1]

# The predicted genes are sequentially named by PRODIGAL in this fashion: contigname_orfnumber
# These two functions retrieve both parts of the name
def contig_name(seqid):
    return '_'.join(seqid.split('_')[:-1])

def orf_number(seqid):
    return seqid.split('_')[-1]

# Three dictionaries are built with a single for loop;
# The first one lists, for each crAss-like contig in the database, the reference sequences which match its own proteins
contigs_hits = {}

# In the second one every query sequence is linked to all the reference sequences it hits
queries_hits = {}

for line in blast_output:
    if contig_name(line[1]) not in contigs_hits:
        contigs_hits[contig_name(line[1])] = {}
    if line[0] not in contigs_hits[contig_name(line[1])]:
        contigs_hits[contig_name(line[1])][line[0]] = []
    if line[1] not in contigs_hits[contig_name(line[1])][line[0]]:
        contigs_hits[contig_name(line[1])][line[0]].append(line[1])

    if line[0] not in queries_hits:
        queries_hits[line[0]] = []
    if line[1] not in queries_hits[line[0]]:
        queries_hits[line[0]].append(line[1])

# The FASTA file containing the predicted genes is parsed
pred_genes = SeqIO.index('/home/alessandro/crassphage_def/all_crassphage_clean.pred.fna', 'fasta')

# crAss-like contigs that had less than two sequences matching with the reference are eliminated
contigs_to_avoid = [contig for contig in contigs_hits if len(contigs_hits[contig]) < 2]

# All the alignment units are created, except for the "special cases" concerning a gene duplication and an insertion
for query in queries_hits:
    if query not in ('NC_024711.1_34', 'NC_024711.1_36', 'NC_024711.1_45', 'NC_024711.1_91'):
        cluster = []
        for hit in queries_hits[query]:
            if contig_name(hit) not in contigs_to_avoid:
                cluster.append(pred_genes[hit])
        SeqIO.write(cluster, 'cluster_' + str(orf_number(query)) + '.fna', 'fasta')

# An individual FASTA file is created with all the sequences matching both parts of the dUTPase protein.
# They will be later split manually
for query in queries_hits:
    if query in ('NC_024711.1_34', 'NC_024711.1_36'):
        cluster = []
        for hit in queries_hits[query]:
            if contig_name(hit) not in contigs_to_avoid:
                cluster.append(pred_genes[hit])
        SeqIO.write(cluster, 'cluster_34_36.fna', 'fasta')

# Three separate alignment units are created for the RepL proteins
# one comprising all the homolog genes
cluster_45_91 = []
for hit in queries_hits['NC_024711.1_45']:
    if contig_name(hit) not in contigs_to_avoid:
        cluster_45_91.append(pred_genes[hit])
SeqIO.write(cluster_45_91, 'cluster_45_91.fna', 'fasta')

# one comprising the sequences that only match with the shorter paralog
cluster_45 = []
for hit in queries_hits['NC_024711.1_45']:
    if contig_name(hit) not in contigs_to_avoid and hit not in queries_hits['NC_024711.1_91']:
        cluster_45.append(pred_genes[hit])
SeqIO.write(cluster_45, 'cluster_45.fna', 'fasta')

# and the last one comprising the sequences matching with the longer paralog
cluster_91 = []
for hit in queries_hits['NC_024711.1_91']:
    if contig_name(hit) not in contigs_to_avoid:
        cluster_91.append(pred_genes[hit])
SeqIO.write(cluster_91, 'cluster_91.fna', 'fasta')
