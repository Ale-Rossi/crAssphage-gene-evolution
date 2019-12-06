"""
A script which averages the information quantity for each alignment and writes all the values to a single file
"""

with open('average_info', 'w') as f:
    for i in range(1, 93):
        filename = 'cluster_' + str(i) + '.csv'
        reference_seq_name = 'NC_024711.1_' + str(i)
        with open(filename, 'r') as g:
            infos = [float(line.strip().split('\t')[0]) for line in g.readlines()[1:]]
        avginf = sum(infos)/len(infos)
        f.write(reference_seq_name + '\t' + str(avginf))
        f.write('\n')
