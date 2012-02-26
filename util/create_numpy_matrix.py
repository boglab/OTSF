import pickle
from math import log
import numpy as np

scoring_matrix_dtype = np.dtype([
    ('HD', 'float32'),
    ('NI', 'float32'),
    ('NG', 'float32'),
    ('NN', 'float32'),
    ('NS', 'float32'),
    ('N*', 'float32'),
    ('HG', 'float32'),
    ('HA', 'float32'),
    ('ND', 'float32'),
    ('NK', 'float32'),
    ('HI', 'float32'),
    ('HN', 'float32'),
    ('NA', 'float32'),
    ('IG', 'float32'),
    ('H*', 'float32'),
    ('S*', 'float32'),
    ('NH', 'float32'),
    ('YG', 'float32'),
    ('SN', 'float32'),
    ('SS', 'float32'),
    ('NC', 'float32'),
    ('HH', 'float32'),
    #other
    ('ZZ', 'float32'),
])

scoring_matrix = np.zeros(4, dtype=scoring_matrix_dtype)

DNA_dict = {'A':0, 'C':1, 'G':2, 'T':3}

comp_DNA_dict = {'A':3, 'C':2, 'G':1, 'T':0}

diresidue_counts = {}

# Diresidue counts for computing entropy
# known
diresidue_counts['HD'] = [	7,	99,	0,	 1]
diresidue_counts['NI'] = [	58,	 6,	0,	 0]
diresidue_counts['NG'] = [	6,	 6,	1,	57]
diresidue_counts['NN'] = [	21,	 8,	26,	 2]
diresidue_counts['NS'] = [	20,	 6,	4,	 0]
diresidue_counts['N*'] = [	1,	11,	1,	 7]
diresidue_counts['HG'] = [	1,	 2,	0,	15]
diresidue_counts['HA'] = [	1,	 4,	1,	 0]
diresidue_counts['ND'] = [	0,	 4,	0,	 0]
diresidue_counts['NK'] = [	0,	 0,	2,	 0]
diresidue_counts['HI'] = [	0,	 1,	0,	 0]
diresidue_counts['HN'] = [	0,	 0,	1,	 0]
diresidue_counts['NA'] = [	0,	 0,	1,	 0]
diresidue_counts['IG'] = [	0,	 0,	0,	 1]
diresidue_counts['H*'] = [	0,	 0,	0,	 1]

# unknown
diresidue_counts['S*'] = [	1,	1,	1,	1]
diresidue_counts['NH'] = [	1,	1,	1,	1]
diresidue_counts['YG'] = [	1,	1,	1,	1]
diresidue_counts['SN'] = [	1,	1,	1,	1]
diresidue_counts['SS'] = [	1,	1,	1,	1]
diresidue_counts['NC'] = [	1,	1,	1,	1]
diresidue_counts['HH'] = [	1,	1,	1,	1]


#known diresidue counts for computing RVD frequencies for each nucleotide
known_diresidues = ['HD', 'NI', 'NG', 'NN', 'NS', 'N*', 'HG', 'HA', 'ND', 'NK', 'HI', 'HN', 'NA', 'IG', 'H*']

# Process diresidue counts based on user-defined weight (default=0.9)

weight = 0.9

diresidue_probability = {}

for diresidue, counts in diresidue_counts.iteritems():
    diresidue_probability[diresidue] = [
        10 * (float(counts[0]) / sum(counts) * weight + (1 - weight) / 4),
        10 * (float(counts[1]) / sum(counts) * weight + (1 - weight) / 4),
        10 * (float(counts[2]) / sum(counts) * weight + (1 - weight) / 4),
        10 * (float(counts[3]) / sum(counts) * weight + (1 - weight) / 4),
    ]

#other
diresidue_probability['ZZ'] = [0.25,	0.25,	0.25,	0.25]


#Normalize probabilities and create scoring matrix
for diresidue, probabilities in diresidue_probability.iteritems():
    total = float(sum(probabilities))
    scoring_matrix[diresidue] = [(-1 * log(probabilities[i] / total)) for i in range(4)]

with open("numpy_matrix_dump", "w") as scoring_matrix_file:
	pickle.dump(scoring_matrix, scoring_matrix_file)
	