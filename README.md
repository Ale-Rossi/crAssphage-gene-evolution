# crAssphage-gene-evolution
A study on the evolution of genes of the crAssphage virus.

Script 01 was used in order to eliminate duplicate contigs from the database.
Script 02 was used to build alignment units from the ORFs predicted by PRODIGAL, according to the BLAST output.
Script 03 was developed to run a pairwise comparison of the Maximum Likelihood matrices output by IQTREE.
Script 04 was used to generate the heatmap representing the Mirrortree coefficient computed iwith the previous script.
Script 05 was used to calculate, for each alignment, the quantity of information for each position.
The quantity of information was then averaged out along the whole sequence length with script 06.
In order to retrieve the annotation for crAssphage ORFs, a BLAST run was executed between the ORFs predicted by PRODIGAL and those used in the study by Yutin and colleagues. The sequences of the two sets were then matched to each other with script 07, and the functional role for each ORF was assigned with script 08.
All the statistics which were computed separatedly were then gathered in a single tabular file with script 09.
Violin plots showing the distribution of various statistics according to functional groups were generated with script 10.
Two more violin plots showing the effect of the filtering steps were produced with scripts 11 and 12.
Script 13 was used to plot the information quantity vs. the number of sequences both in the crAss-like alignment unit produced in this study and the protein viral families from the pVOGs database.
