"""
This script assigns a functional module to each gene, based on the annotation retrieved from the study by Yutin et al., 2018
"""

with open('prodigal_yutin_match', 'r') as f:
    tab = [line.strip().split('\t') for line in f.readlines()]

for i in range(len(tab)):
    if tab[i][2] == 'Uncharacterized':
        tab[i].append('Uncharacterized')
    elif i <= 42:
        tab[i].append('Replication')
    elif i <= 44:
        tab[i].append('Transcription')
    elif i <= 72:
        tab[i].append('Tail and structural protein')
    elif i <= 76:
        tab[i].append('Capsid')
    elif i <= 86:
        tab[i].append('Other')
    else:
        tab[i].append('Replication')

new_tab = []
for i in range(1, 93):
    orf = 'NC_024711.1_' + str(i)
    if orf in [line[0] for line in tab]:
        for line in tab:
            if line[0] == orf:
                new_tab.append(line)
    else:
        new_tab.append([orf, 'none', 'Uncharacterized', 'Uncharacterized'])

with open('prodigal_yutin_match_II', 'w') as f:
    for line in new_tab:
        f.write('\t'.join(line) + '\n')
