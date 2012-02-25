from base_task import TalentTask, _workerState
from celery.signals import worker_init
from genome_scoring import GENOME_FILE, TREE_FILE, BOUNDARY_FILE
from sstree import PySSTree
from blist import sortedlist
import pickle

def _initializeWorker(**kwargs):
    
    print '_initializeWorker called'
    
    organism = "mus_musculus"
    
    with open(GENOME_FILE % organism, 'r') as f:
        text = f.read()

    _workerState["sTree"] = PySSTree(text, "load", TREE_FILE % organism)
    
    with open(BOUNDARY_FILE % organism, "r") as f:
        _workerState["geneBoundaries"] = sortedlist(pickle.load(f), key=lambda x: x['pos'])
    
worker_init.connect(_initializeWorker)
