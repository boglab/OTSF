import pickle
from math import log
import numpy as np

# X replaces * in RVD names
# The scorer replaces any * in the RVD sequence with X
# Z would have been used instead of X, but someone defined HZ as a macro somewhere in the cython compiler software stack

scoring_matrix = np.zeros((23, 4), dtype=np.float32, order='F')

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
diresidue_counts['NX'] = [	1,	11,	1,	 7]
diresidue_counts['HG'] = [	1,	 2,	0,	15]
diresidue_counts['HA'] = [	1,	 4,	1,	 0]
diresidue_counts['ND'] = [	0,	 4,	0,	 0]
diresidue_counts['NK'] = [	0,	 0,	2,	 0]
diresidue_counts['HI'] = [	0,	 1,	0,	 0]
diresidue_counts['HN'] = [	0,	 0,	1,	 0]
diresidue_counts['NA'] = [	0,	 0,	1,	 0]
diresidue_counts['IG'] = [	0,	 0,	0,	 1]
diresidue_counts['HX'] = [	0,	 0,	0,	 1]

# unknown
diresidue_counts['SX'] = [	1,	1,	1,	1]
diresidue_counts['NH'] = [	1,	1,	1,	1]
diresidue_counts['YG'] = [	1,	1,	1,	1]
diresidue_counts['SN'] = [	1,	1,	1,	1]
diresidue_counts['SS'] = [	1,	1,	1,	1]
diresidue_counts['NC'] = [	1,	1,	1,	1]
diresidue_counts['HH'] = [	1,	1,	1,	1]

numpyMap = {}

numpyMap['HD'] = 0
numpyMap['NI'] = 1
numpyMap['NG'] = 2
numpyMap['NN'] = 3
numpyMap['NS'] = 4
numpyMap['NX'] = 5
numpyMap['HG'] = 6
numpyMap['HA'] = 7
numpyMap['ND'] = 8
numpyMap['NK'] = 9
numpyMap['HI'] = 10
numpyMap['HN'] = 11
numpyMap['NA'] = 12
numpyMap['IG'] = 13
numpyMap['HX'] = 14
numpyMap['SX'] = 15
numpyMap['NH'] = 16
numpyMap['YG'] = 17
numpyMap['SN'] = 18
numpyMap['SS'] = 19
numpyMap['NC'] = 20
numpyMap['HH'] = 21
numpyMap['XX'] = 22


#known diresidue counts for computing RVD frequencies for each nucleotide
known_diresidues = ['HD', 'NI', 'NG', 'NN', 'NS', 'NX', 'HG', 'HA', 'ND', 'NK', 'HI', 'HN', 'NA', 'IG', 'HX']

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
diresidue_probability['XX'] = [0.25,	0.25,	0.25,	0.25]


#Normalize probabilities and create scoring matrix
for diresidue, probabilities in diresidue_probability.iteritems():
    total = float(sum(probabilities))
    numpyIndex = numpyMap[diresidue]
    scoring_matrix[numpyIndex] = [(-1 * log(probabilities[i] / total)) for i in range(4)]

with open("numpy_scoring_dump_2d", "w") as scoring_matrix_file:
	pickle.dump(scoring_matrix, scoring_matrix_file)
	