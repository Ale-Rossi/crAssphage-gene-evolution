"""
This script was used to retrieve all the data about each crAssphage gene and compile them in a single tabular file
"""

orf_names = ['NC_024711.1_' + str(i) for i in range(1,93)]

annotations = {}
category = {}
with open('prodigal_yutin_match_II', 'r') as f:
    for line in f.readlines():
        line = line.strip().split('\t')
        if line[0] in orf_names:
            annotations[line[0]] = line[2]
            category[line[0]] = line[3]

num_seqs_in_alignment = {}
with open('/home/alessandro/controllare_per_la_tesi/aligned_2.0/num_seqs_in_alignments.txt', 'r') as f:
    for line in f.readlines():
        line = line.strip().split('\t')
        if line[0] in orf_names:
            num_seqs_in_alignment[line[0]] = line[1]

avg_info_quant = {}
with open('/home/alessandro/controllare_per_la_tesi/info_quantity_profiles/average_info', 'r') as f:
    for line in f.readlines():
        line = line.strip().split('\t')
        if line[0] in orf_names:
            avg_info_quant[line[0]] = line[1]

avg_pearson = {}
with open('/home/alessandro/controllare_per_la_tesi/matrices/pearson_avg', 'r') as f:
    for line in f.readlines():
        line = line.strip().split('\t')
        if line[0] in orf_names:
            avg_pearson[line[0]] = line[1]

with open('/home/alessandro/controllare_per_la_tesi/all_comprehensive_tab.csv', 'w') as f:
    f.write('\t'.join(['NAME', 'ANNOTATION', 'CATEGORY', 'NR. OF SEQUENCES', 'AVERAGE INFORMATION QUANTITY', 'AVERAGE PEARSON CORRELATION']))
    f.write('\n')
    for orf in orf_names:
        f.write('\t'.join([orf, annotations[orf], category[orf], num_seqs_in_alignment[orf], avg_info_quant[orf], avg_pearson[orf]]))
        f.write('\n')
