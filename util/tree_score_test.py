import pickle

from genome_scoring import SCORING_MATRIX, GENOME_FILE, TREE_FILE, BOUNDARY_FILE
from sstree import PySSTree, ScoreTalentTask
from blist import sortedlist

import string
import numpy as np

import time
import timeit

import random

ORGANISM = "arabidopsis_thaliana"

scoringMatrix = pickle.load(open(SCORING_MATRIX, "r"))

with open(GENOME_FILE % ORGANISM, 'r') as f:
        text = f.read()

sTree = PySSTree(text, "load", TREE_FILE % ORGANISM)

with open(BOUNDARY_FILE % ORGANISM, "r") as f:
    geneBoundaries = sortedlist(pickle.load(f), key=lambda x: x['pos'])
    
subMatrix = {
    "NG": ["T"],
    "NN": ["A", "G"],
    "HD": ["C"],
    "NI": ["A"],
    "NS": ["A"]
}


baseMap = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}


revBaseMap = {
    'A': 3,
    'C': 2,
    'G': 1,
    'T': 0
}


randMap = ["NG", "NN", "HD", "NI", "NS"]

#querySequence = [randMap[random.randrange(0, 5, 1)] for i in range(1,22)]
querySequence = "NN NG NN NG NG HD HD NG HD NI NN HD NG NN NG NN NG".split()

def doIt():
        ScoreTalentTask(string.join(querySequence, " "), "droso_standalone_seq_" + str(int(time.time())) + ".txt", subMatrix, baseMap, geneBoundaries, scoringMatrix, sTree)
        
print timeit.Timer(doIt).timeit(1)