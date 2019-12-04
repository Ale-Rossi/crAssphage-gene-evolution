
# coding: utf-8

# In[1]:


with open('/home/alessandro/crassphage_def/all_crassphage_clean.blast') as f:
    blast_output = [line.strip().split('\t') for line in f.readlines() if len(line.strip().split('\t')) > 1]


# In[2]:


def contig_name(seqid):
    return '_'.join(seqid.split('_')[:-1])


# In[3]:


def orf_number(seqid):
    return seqid.split('_')[-1]


# In[4]:


contigs_hits = {}

for line in blast_output:
    if contig_name(line[1]) not in contigs_hits:
        contigs_hits[contig_name(line[1])] = {}
    if line[0] not in contigs_hits[contig_name(line[1])]:
        contigs_hits[contig_name(line[1])][line[0]] = []
    if line[1] not in contigs_hits[contig_name(line[1])][line[0]]:
        contigs_hits[contig_name(line[1])][line[0]].append(line[1])


# In[5]:


queries_hits = {}

for line in blast_output:
    if line[0] not in queries_hits:
        queries_hits[line[0]] = []
    if line[1] not in queries_hits[line[0]]:
        queries_hits[line[0]].append(line[1])


# In[6]:


hits_queries = {}

for line in blast_output:
    if line[1] not in hits_queries:
        hits_queries[line[1]] = []
    if line[0] not in hits_queries[line[1]]:
        hits_queries[line[1]].append(line[0])


# In[7]:


from Bio import SeqIO


# In[8]:


pred_genes = SeqIO.index('/home/alessandro/crassphage_def/all_crassphage_clean.pred.fna', 'fasta')


# In[9]:


contigs_to_avoid = [contig for contig in contigs_hits if len(contigs_hits[contig]) < 2]


# In[10]:


for query in queries_hits:
    if query not in ('NC_024711.1_34', 'NC_024711.1_36', 'NC_024711.1_45', 'NC_024711.1_91'):
        cluster = []
        for hit in queries_hits[query]:
            if contig_name(hit) not in contigs_to_avoid:
                cluster.append(pred_genes[hit])
        SeqIO.write(cluster, 'cluster_' + str(orf_number(query)) + '.fna', 'fasta')


# In[11]:


for query in queries_hits:
    if query in ('NC_024711.1_34', 'NC_024711.1_36'):
        cluster = []
        for hit in queries_hits[query]:
            if contig_name(hit) not in contigs_to_avoid:
                cluster.append(pred_genes[hit])
        SeqIO.write(cluster, 'cluster_34_36.fna', 'fasta')


# In[12]:


cluster_45_91 = []
for hit in queries_hits['NC_024711.1_45']:
    if contig_name(hit) not in contigs_to_avoid:
        cluster_45_91.append(pred_genes[hit])
SeqIO.write(cluster_45_91, 'cluster_45_91.fna', 'fasta')


# In[13]:


cluster_45 = []
for hit in queries_hits['NC_024711.1_45']:
    if contig_name(hit) not in contigs_to_avoid and hit not in queries_hits['NC_024711.1_91']:
        cluster_45.append(pred_genes[hit])
SeqIO.write(cluster_45, 'cluster_45.fna', 'fasta')


# In[14]:


cluster_91 = []
for hit in queries_hits['NC_024711.1_91']:
    if contig_name(hit) not in contigs_to_avoid:
        cluster_91.append(pred_genes[hit])
SeqIO.write(cluster_91, 'cluster_91.fna', 'fasta')

