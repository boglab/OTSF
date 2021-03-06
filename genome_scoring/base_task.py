from celery.task import task
from celery.registry import tasks
from celery.signals import worker_init

import pickle
import re
import numpy as np

from sstree import ScoreTalentTask
from genome_scoring import SCORING_MATRIX, GENOME_FILE, TREE_FILE, BOUNDARY_FILE, DRUPAL_CALLBACK_URL

import urllib

_workerState = {}

_workerState["scoringMatrix"] = np.asfortranarray(pickle.load(open(SCORING_MATRIX, "r")))

@task
def TalentTask(querySequence, outputFilepath, nodeID, revComp):
    ScoreTalentTask(querySequence, outputFilepath, revComp, _workerState["geneBoundaries"], _workerState["scoringMatrix"], _workerState["sTree"])
    if nodeID != -1:
        urllib.urlopen(DRUPAL_CALLBACK_URL + str(nodeID) + "/0")

