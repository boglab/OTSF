from celery.task import task
from celery.registry import tasks
from celery.signals import worker_init

import pickle
import re

from sstree import ScoreTalentTask
from genome_scoring import SCORING_MATRIX, GENOME_FILE, TREE_FILE, BOUNDARY_FILE, DRUPAL_CALLBACK_URL

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
def TalentTask(querySequence, outputFilepath, nodeID, revComp):#, subMatrix, baseMap, scoringMatrix, sTree):
    ScoreTalentTask(querySequence, outputFilepath, revComp, _workerState["geneBoundaries"], _workerState["scoringMatrix"], _workerState["sTree"])
    if nodeID != -1:
        urllib.urlopen(DRUPAL_CALLBACK_URL + str(nodeID) + "/0")

