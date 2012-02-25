from celery.task import task
from celery.registry import tasks
from celery.signals import worker_init

import pickle
import re

from sstree import ScoreTalentTask
from genome_scoring import SCORING_MATRIX, GENOME_FILE, TREE_FILE, BOUNDARY_FILE

import urllib

_workerState = {}

_workerState["subMatrix"] = {
    "NG": ["T"],
    "NN": ["A", "G"],
    "HD": ["C"],
    "NI": ["A"],
    "NS": ["A"]
};

_workerState["baseMap"] = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
};

_workerState["revBaseMap"] = {
    'A': 3,
    'C': 2,
    'G': 1,
    'T': 0
};

_workerState["scoringMatrix"] = pickle.load(open(SCORING_MATRIX, "r"))


@task
def TalentTask(querySequence, outputFilepath, nodeID):#, subMatrix, baseMap, scoringMatrix, sTree):
    ScoreTalentTask(querySequence, outputFilepath, _workerState["subMatrix"], _workerState["baseMap"], _workerState["geneBoundaries"], _workerState["scoringMatrix"], _workerState["sTree"])
    if nodeID != -1:
        urllib.urlopen("https://boglab.plp.iastate.edu/talent/jobcomplete/" + str(nodeID) + "/0")

