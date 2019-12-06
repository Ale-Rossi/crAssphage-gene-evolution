"""
A script for calculating the Mirrortree coefficient for each pair of crAssphage genes.
It compares the Maximum Likelihood matrices by trimming them, reordering them and converting them to vectors, then computing the Pearson's correlation coefficient between the pair of vectors
It is a very clunky script which could be replaced by already existing software, such as the Mantel test found in Scikit-bio,
or just benefit from the implementation of concurrency (most easily by running it with an external program, such as GNU Parallel)
"""

from os import listdir
from numpy import mat as numpy_mat
from scipy.spatial.distance import squareform, correlation, euclidean
from scipy.stats import pearsonr, spearmanr
from copy import deepcopy
from itertools import combinations

# Each matrix file is read into a list of lists
matrices_names = [filename for filename in listdir() if 'mldist' in filename]
matrices = {}
for filename in matrices_names:
    with open(filename, 'r') as file:
        matrices[filename] = [line.strip().split(' ') for line in file.readlines()[1:]] # As the first line of each matrix only contains the number of leaves in the tree, it is ignored 

# Sequence names are contig names followed by an underscore and the number of the ORF. This function returns the contig name of sequences.
def contigname(seqname):
    return '_'.join(seqname.split('_')[:-1])

# This function substitutes the name of the sequence with the name of the relative contig, in order for the matrices to be compared
def givecontigs(matrix):
    for line in matrix:
        line[0] = contigname(line[0])

# The matrices contain many empty strings. This function removes them.
def cleanmatrix(matrix):
    for line in matrix:
        while '' in line:
            line.remove('')

# This function turns every element of each row of a matrix into a tuple consisting of the name of the second term of comparison and the value itself
# This allows for each row to be reordered alphabetically
def annotate_matrix(matrix):
    indexes = [row[0] for row in matrix]
    for row in matrix:
        for i in range(len(indexes)):
            row[i+1] = (indexes[i], row[i+1])

# Finally the matrix is reordered.
def reorder_matrix(matrix):
    matrix.sort(key= lambda row: row[0])		# Rows are sorted according to their first element
    for i in range(len(matrix)):
        matrix[i] = [matrix[i][0], matrix[i][1:]]	# In each row the values are separated from the index, in order for them to be sorted separately form the index
        for row in matrix:				# The values of each row are sorted using the first element of the tuple
        row[1].sort(key= lambda element : element[0])

# The previously defined functions are thus enacted on each matrix
for matrix in matrices:
    givecontigs(matrices[matrix])
    cleanmatrix(matrices[matrix])
    annotate_matrix(matrices[matrix])
    reorder_matrix(matrices[matrix])

# This function is executed when comparing two matrices.
# It trims a matrix row-wise, keeping only the rows found in a specific list
def del_mat_rows(matrix, namelist):
    return [row for row in matrix if row[0] in namelist]

# Similarly, this function trims each row element-wise.
def del_mat_cols(matrix, namelist):
    matrix_out = deepcopy(matrix)
    for row in matrix_out:
        row[1] = [tuple for tuple in row[1] if tuple[0] in namelist]
    return matrix_out

# Another function restores the matrix, leaving only the numeric values
def matrix_restorer(matrix):
    matrix_out = deepcopy(matrix)
    for row in matrix_out:
        for i in range(len(row[1])):
            row[1][i] = float(row[1][i][1])
    for i in range(len(matrix_out)):
        matrix_out[i] = matrix_out[i][1]
    return matrix_out

# The pairwise comparisons are carried out and the results are written on a .csv file
# Other coefficients are calculated too

with open('matrix_correlations.csv', 'w') as output_file:
    
    output_file.write('\t'.join(['matr_1', 'matr_2', 'nr. shared contigs', 'spearman correlation', 'pearson correlation', 'euclidean distance', 'correlation distance']))
    output_file.write('\n')

    for matrix_name_a, matrix_name_b in combinations(matrices_names, 2):
    
        mat_a = deepcopy(matrices[matrix_name_a])
        mat_b = deepcopy(matrices[matrix_name_b])
        print(matrix_name_a, matrix_name_b)
        
        names_a = [row[0] for row in mat_a]
        names_b = [row[0] for row in mat_b]
        common_names = [x for x in names_a if x in names_b]
        
        if len(common_names) < 2:				# If the matrix pair has less than two elements in common, it is not possible to compare them.
            
            output_file.write('\t'.join([matrix_name_a, matrix_name_b, str(len(common_names)), 'NA', 'NA', 'NA', 'NA']))
            output_file.write('\n')
            
        else:
            mat_a = del_mat_rows(mat_a, common_names)
            mat_a = del_mat_cols(mat_a, common_names)
            mat_a = matrix_restorer(mat_a)
            mat_b = del_mat_rows(mat_b, common_names)
            mat_b = del_mat_cols(mat_b, common_names)
            mat_b = matrix_restorer(mat_b)           
            symmat_a = numpy_mat(mat_a)				# The matrices are coverted into numpy matrices
            symmat_b = numpy_mat(mat_b)
            vec_a = squareform(symmat_a, checks=False)
            vec_b = squareform(symmat_b, checks=False)
            corr_dist = correlation(vec_a, vec_b)
            pears = pearsonr(vec_a, vec_b)
            spear = spearmanr(vec_a, vec_b, axis= None)
            eucl = euclidean(vec_a, vec_b)
            print(spear)
            print(pears)
            print(eucl)
            print(corr_dist)
            output_file.write('\t'.join([matrix_name_a, matrix_name_b, str(len(common_names)), str(spear), str(pears), str(eucl), str(corr_dist)]))
            output_file.write('\n')
