"""
This script was developed in order to create a scatter plot showing the relation between the abundance of a certain proteic family and the sequence conservation, represented by the average quantity of information of the alignment.
It makes use of the pyplot package
"""

import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt

podo_vogs = pd.read_csv("average_info_pop_pvogs", sep='\t')

vogs_sizes = pd.read_csv("vogs_sizes", sep='\t')

podo_vogs_complete = pd.merge_ordered(podo_vogs, vogs_sizes, on='VOG')

crass_like = pd.read_csv('../controllare_per_la_tesi/work_with_statistics/all_comprehensive_tab.csv', sep='\t')

plt.rcParams['figure.figsize'] = 10, 10
plt.plot(podo_vogs_complete['size'], podo_vogs_complete['average_info'], 'o', crass_like['NR. OF SEQUENCES'], crass_like['AVERAGE INFORMATION QUANTITY'], 'ro')
plt.xscale('log')
plt.xlabel("Number of sequences")
plt.ylabel("Quantity of information (bits)")
plt.legend(("pVOGs", "crAss-like alignment units"))

log_size_podo = pd.Series([np.log10(x) for x in podo_vogs_complete['size']])
log_size_crass = pd.Series([np.log10(x) for x in crass_like['NR. OF SEQUENCES']])

stats.linregress(log_size_podo, podo_vogs_complete['average_info'])
stats.linregress(log_size_crass, crass_like['AVERAGE INFORMATION QUANTITY'])
