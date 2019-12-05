"""
This script computes the average Pearson's correlation coefficient for each phylogenetic tree.
It also makes use of the pyplot framework for generating two of the plots presented in the paper
"""

from matplotlib import pyplot as plt
from os import listdir

matrices_names = [filename for filename in listdir() if 'mldist' in filename]

# Each Maximum Likelihood matrix is related to the phylogenetic tree of a specific crAssphage ORF.
# Each ORF is sequentially numbered.
# The following function retrieves the ORF number from the name of the matrix file

def matrix_num(matrix_name):
    if 'clean' in matrix_name:
        return int(matrix_name.split('_')[1])
    elif 'fasta' in matrix_name:
        return int(matrix_name.replace('cluster_', '').replace('.fasta.mldist', ''))
    else:
        return int(matrix_name.replace('cluster_', '').replace('.msa.mldist', ''))

matrices_names.sort(key = matrix_num)

# The .csv file containing the pairwise correlation coefficients is read
with open('matrix_correlations.csv', 'r') as f:
    correlations = [line.strip().split('\t') for line in f.readlines()[1:]] # The header is skipped

# The Pearson's correlation coefficient is written as (coefficient, p-value)
# This function extracts the coefficient
def extract_pearson_corr_coef(x):
    return float(x.split(',')[0].replace('(', ''))

# The square matrix containing the Mirrortree coefficient of each pair of phylogenetic trees is created
pearson_square = []
for i in range(len(matrices_names)):
    pearson_square.append([])
    
for line in pearson_square:
    for i in range(len(matrices_names)):
        line.append(-100)

# -100 is given as a default value for each pair of matrices
# It is replaced by the actual correlation coefficient except when the trees have less than 5 sequences in common
# So that they can easily be painted white in the plot

for line in correlations:
    if int(line[2]) > 5:
        i = matrices_names.index(line[0])
        j = matrices_names.index(line[1])
        pearson_square[i][j] = extract_pearson_corr_coef(line[4])
        pearson_square[j][i] = extract_pearson_corr_coef(line[4])

# The square matrix is also written to a .csv file

with open('pearson_square.csv', 'w') as f:
    for line in pearson_square:
        f.write('\t'.join([str(x) for x in line]) + '\n')

# The heatmap is generated and saved. The colormap is customized in order to paint coefficients < -1 (i.e. those between trees which shared less than 5 sequences) in white

plt.rcParams['figure.figsize'] = (12, 10)
custom_cmap = plt.get_cmap('plasma')
custom_cmap.set_under('w')
plt.pcolormesh(pearson_square, vmin = -1, vmax = 1, cmap = custom_cmap)
plt.colorbar()
plt.savefig('heatmap_pearson.jpg')

# All the correlations are saved in a list, which is then used to generate the histogram

pearson_correlations = []
for line in correlations:
    if line[0] != line[1] and int(line[2]) > 5:
        pearson_correlations.append(extract_pearson_corr_coef(line[4]))

plt.rcParams['figure.figsize'] = 12, 7
plt.hist(pearson_correlations, bins=50)
plt.xlabel('pearson correlation coefficient')
plt.ylabel('Number of instances')
plt.savefig('pearson_correlation_histogram.jpg')

# The average Mirrortree coefficient for each tree is written in another file

with open('pearson_avg', 'w') as f:
    for i in range(len(matrices_names)):
        f.write('NC_024711.1_' + str(matrix_num(matrices_names[i])) + '\t' + str(average_pearson_correlation[i]))
        f.write('\n')
