import pickle
from optparse import OptionParser
import sys


DNA_dict = {'A':0, 'C':1, 'G':2, 'T':3}
subMatrix = {
    "NG": ["T"],
    "NN": ["A", "G"],
    "HD": ["C"],
    "NI": ["A"],
    "NS": ["A"]
};


scoringMatrix = pickle.load(open("scoring_matrix_dump", "r"))


dnaSequence = sys.argv[1]
rvdSequence = sys.argv[2]

score = 0
mismatches = 0

for i,dr in enumerate(rvdSequence.split()):
    base = dnaSequence[i]
    matrixPos = DNA_dict[base]
    if base not in subMatrix[dr]:
        mismatches += 1
    score+= scoringMatrix[dr][matrixPos]

print "mismatches: %d" % (mismatches)
print "score: %lf" % (score)
