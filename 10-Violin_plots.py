"""
This script was developed in order to visualize data according to functional categories with the pyplot package
It makes use of the data previously gathered in the tabular file "all_comprehensive_tab.csv"
"""

import matplotlib.pyplot as plt
from scipy import stats

names = []
annotations = []
categories = []
cluster_dimensions = []
avg_info_quant = []
avg_spearman = []
avg_pearson = []
labelpos = range(1,7)
labels = ('Uncharacterized', 'Replication', 'Transcription', 'Tail and structural', 'Capsid', 'Other')


with open('all_comprehensive_tab.csv', 'r') as f:
    tab = [line.strip().split('\t') for line in f.readlines()[1:]]
    for line in tab:
        names.append(line[0])
        annotations.append(line[1])
        categories.append(line[2])
        cluster_dimensions.append(int(line[3]))
        avg_info_quant.append(float(line[4]))
        avg_spearman.append(float(line[5]))
        avg_pearson.append(float(line[6]))


# The first plot represents the average information according to functional category

uncharacterized_info = []
replication_info = []
trascription_info = []
tail_struc_info = []
capsid_info = []
other_info = []
average_infos_according_category = [uncharacterized_info, replication_info, trascription_info, tail_struc_info, capsid_info, other_info]
 
for i in range(len(categories)):
    if categories[i] == 'Uncharacterized':
        uncharacterized_info.append(avg_info_quant[i])
    elif categories[i] == 'Replication':
        replication_info.append(avg_info_quant[i])
    elif categories[i] == 'Transcription':
        trascription_info.append(avg_info_quant[i])
    elif categories[i] == 'Tail and structural protein':
        tail_struc_info.append(avg_info_quant[i])
    elif categories[i] == 'Capsid':
        capsid_info.append(avg_info_quant[i])
    elif categories[i] == 'Other':
        other_info.append(avg_info_quant[i])

plt.violinplot(average_infos_according_category, showmedians=True)
plt.xticks(labelpos, labels)
plt.ylabel('Average information quantity')


# The second plot represents the average Mirrortree coefficient according to functional category

uncharacterized_pearson = []
replication_pearson = []
trascription_pearson = []
tail_struc_pearson = []
capsid_pearson = []
other_pearson = []
average_pearson_according_category = [uncharacterized_pearson, replication_pearson, trascription_pearson, tail_struc_pearson, capsid_pearson, other_pearson]
 
for i in range(len(categories)):
    if categories[i] == 'Uncharacterized':
        uncharacterized_pearson.append(avg_pearson[i])
    elif categories[i] == 'Replication':
        replication_pearson.append(avg_pearson[i])
    elif categories[i] == 'Transcription':
        trascription_pearson.append(avg_pearson[i])
    elif categories[i] == 'Tail and structural protein':
        tail_struc_pearson.append(avg_pearson[i])
    elif categories[i] == 'Capsid':
        capsid_pearson.append(avg_pearson[i])
    elif categories[i] == 'Other':
        other_pearson.append(avg_pearson[i])

plt.violinplot(average_pearson_according_category, showmedians = True)
plt.xticks(labelpos, labels)
plt.ylabel('Average pearson correlation coefficient')


# In the third plot the quantity of sequences in each alignment is represented

uncharacterized_clusdim = []
replication_clusdim = []
trascription_clusdim = []
tail_struc_clusdim = []
capsid_clusdim = []
other_clusdim = []
clusdim_according_category = [uncharacterized_clusdim, replication_clusdim, trascription_clusdim, tail_struc_clusdim, capsid_clusdim, other_clusdim]
for i in range(len(categories)):
    if categories[i] == 'Uncharacterized':
        uncharacterized_clusdim.append(cluster_dimensions[i])
    elif categories[i] == 'Replication':
        replication_clusdim.append(cluster_dimensions[i])
    elif categories[i] == 'Transcription':
        trascription_clusdim.append(cluster_dimensions[i])
    elif categories[i] == 'Tail and structural protein':
        tail_struc_clusdim.append(cluster_dimensions[i])
    elif categories[i] == 'Capsid':
        capsid_clusdim.append(cluster_dimensions[i])
    elif categories[i] == 'Other':
        other_clusdim.append(cluster_dimensions[i])

plt.violinplot(clusdim_according_category, showmedians=True)
plt.ylabel('Number of sequences in alignment unit')
plt.xticks(labelpos, labels)


# A histogram representing the distribution of the average correlation coefficient for each tree is also created

plt.hist(avg_pearson, bins=30)
plt.xlabel("Average Pearson's correlation coefficient")
plt.ylabel("Nr. of instances")
