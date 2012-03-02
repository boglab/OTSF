import pickle
from optparse import OptionParser
import sys


DNA_dict = {'A':0, 'C':1, 'G':2, 'T':3}
c_DNA_dict = {'A':3, 'C':2, 'G':1, 'T':0}

subMatrix = {
    "NG": ["T"],
    "NN": ["A", "G"],
    "HD": ["C"],
    "NI": ["A"],
    "NS": ["A"]
};


scoringMatrix = pickle.load(open("../genome_scoring/scoring_matrix_dump", "r"))


def reverseComplement(sequence):
    new_sequence = ""
    
    for i in range(len(sequence)):
        base = sequence[len(sequence) - i - 1]
        if base == 'A':
            new_sequence += 'T'
        elif base == 'C':
            new_sequence += 'G'
        elif base == 'G':
            new_sequence += 'C'
        elif base == 'T':
            new_sequence += 'A'
    
    return new_sequence

revcomp = False



if len(sys.argv) == 4 and sys.argv[1] == '-r':
    revcomp = True
    dnaSequence = reverseComplement(sys.argv[2])
    rvdSequence = sys.argv[3]
    diresidues = rvdSequence.split()
    diresidues.reverse()
else:
    dnaSequence = sys.argv[1]
    rvdSequence = sys.argv[2]
    diresidues = rvdSequence.split()

score = 0
mismatches = 0

for i,dr in enumerate(diresidues):
    base = dnaSequence[i]
    if revcomp:
        matrixPos = c_DNA_dict[base]
    else:
        matrixPos = DNA_dict[base]
        
    if base not in subMatrix[dr]:
        mismatches += 1
    score+= scoringMatrix[dr][matrixPos]

print "mismatches: %d" % (mismatches)
print "score: %lf" % (score)
