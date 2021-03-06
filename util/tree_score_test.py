import pickle

#from os.path import abspath, dirname
#import sys
#sys.path.append(dirname(dirname(abspath(__file__))))

from genome_scoring import SCORING_MATRIX, GENOME_FILE, TREE_FILE, BOUNDARY_FILE
from sstree import PySSTree, ScoreTalentTask
from blist import sortedlist

import string
import numpy as np

import time
import timeit

#import random

#ORGANISM = "arabidopsis_thaliana"
ORGANISM = "drosophila_melanogaster"

scoringMatrix = np.asfortranarray(pickle.load(open(SCORING_MATRIX, "r")))

with open(GENOME_FILE % ORGANISM, 'r') as f:
        text = f.read()

sTree = PySSTree(text, "load", TREE_FILE % ORGANISM)

with open(BOUNDARY_FILE % ORGANISM, "r") as f:
    geneBoundaries = sortedlist(pickle.load(f), key=lambda x: x['pos'])

#randMap = ["NG", "NN", "HD", "NI", "NS"]

#querySequence = [randMap[random.randrange(0, 5, 1)] for i in range(1,22)]
#querySequence = "NN NG NN NG NG HD HD NG HD NI NN HD NG NN NG NN NG".split()
querySequence = "HD NG NN NN NI NN NN NG NG HD HD HD NI".split()

def doIt():
        ScoreTalentTask(string.join(querySequence, " "), "test_score_output/droso_standalone_seq_" + str(int(time.time())), True, geneBoundaries, scoringMatrix, sTree)
print "Files loaded, starting test"
print timeit.Timer(doIt).timeit(1)