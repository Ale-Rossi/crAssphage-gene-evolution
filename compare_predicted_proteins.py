"""
The annotation for crAssphage proteins is retrieved from a previous study: Yutin et al., 2018
A BLAST search has been previously been run with the PRODIGAL-predicted ORFs as queries and the ORFs used in the previous study
The ORFs were considered matching if the sequence identity was higher than 90%
"""

with open('prodigal_vs_yutin', 'r') as f:
    blastp = [line.strip().split('\t') for line in f.readlines()]

with open('reference_yutin_annotation.csv') as f:
    yut_ann = [line.strip().split('\t') for line in f.readlines()[1:]]

with open('prodigal_yutin_match', 'w') as f:
    for line in blastp:
        if float(line[2]) > 90:
            f.write('\t'.join([line[0], line[1]]))
            for linee in yut_ann:
                if line[1] in linee:
                    f.write('\t' + linee[2])
            f.write('\n')
