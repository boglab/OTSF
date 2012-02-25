from celery.task import task
from celery.registry import tasks

import pickle
import re

from math import floor
from collections import deque
from bitstring import BitArray

from genome_scoring import SCORING_MATRIX

import string
import urllib

_workerState = {
    "subMatrix": {
        "NG": ["T"],
        "NN": ["A", "G"],
        "HD": ["C"],
        "NI": ["A"],
        "NS": ["A"]
    },
    "baseMap": {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    },
    "revBaseMap": {
        'A': 3,
        'C': 2,
        'G': 1,
        'T': 0
    },
    "scoringMatrix": pickle.load(open(SCORING_MATRIX, "r"))
}



@task
def TalentTask(querySequence, outputFilepath, nodeID):#, subMatrix, baseMap, scoringMatrix, sTree):
    
    scoringMatrix = _workerState["scoringMatrix"]
    subMatrix = _workerState["subMatrix"]
    sTree = _workerState["sTree"]
    baseMap = _workerState["baseMap"]
    revBaseMap = _workerState["revBaseMap"]
    geneBoundaries = _workerState["geneBoundaries"]
    
    querySequence = querySequence.upper().rstrip()
    diresidues = querySequence.split(' ')
    
    bestScore = 0
    for diresidue in diresidues:
        maxProbIndex = min(range(4), key=lambda index: scoringMatrix[diresidue][index])
        bestScore += scoringMatrix[diresidue][maxProbIndex]
    
    cutoffScore = 2.5 * bestScore
    allowedMismatches = floor(0.5 * len(diresidues))
    
    DNARe = re.compile(r'[^ATCG]')
    
    startNode = sTree.search("T")
    
    openSet = deque([startNode])
    scores = {startNode : 0}
    mismatches = {startNode : 0}
    
    with open(outputFilepath, "w") as outputFile:
        
        tableIgnores = []
        #if not revComp:
        tableIgnores.append("Plus strand sequence")
        if len(tableIgnores) > 0:
            outputFile.write("table_ignores:" + string.join(tableIgnores, ",") + "\n")
        outputFile.write("Best Possible Score:" + str(round(bestScore, 2)) + "\n")
        outputFile.write('Genome Coordinates\tStrand\tScore\tTarget Sequence\tPlus strand sequence\n')
        while len(openSet) > 0:
            
            node = openSet.pop()
            parentDepth = sTree.depth(node) - 1
                
            for child in sTree.children(node):
                
                scores[child] = scores[node]
                mismatches[child] = mismatches[node]
                
                if parentDepth < len(diresidues):
                    
                    k = 1
                    depth = parentDepth
                    edgeChar = sTree.edgeChar(child, k)
                    
                    badChar = False
                    
                    while depth < len(diresidues) and edgeChar != '\x00':
                        
                        if DNARe.match(edgeChar):
                            badChar = True
                            break
                        
                        diresidue = diresidues[depth]
                        baseIndex = baseMap[edgeChar]
                        scores[child] += scoringMatrix[diresidue][baseIndex]
                        
                        if edgeChar not in subMatrix[diresidue]:
                            mismatches[child] += 1
                            
                        k += 1
                        depth += 1
                        edgeChar = sTree.edgeChar(child, k)
                    
                    if not badChar and scores[child] < cutoffScore and mismatches[child] <= allowedMismatches:
                        openSet.append(child)
                    
                else:
                    
                    if not sTree.isleaf(child):
                        openSet.append(child)
                    
                    else:
                        
                        textPos = sTree.textpos(child) + 1
                        listIndex = geneBoundaries.bisect_right({'pos': textPos})
                        
                        if listIndex != len(geneBoundaries):
                            
                            if listIndex == 0:
                                outputTextPos = str(textPos)
                            else:
                                outputTextPos = str(textPos - geneBoundaries[listIndex - 1]['pos'] - 1)
                            
                            sequence = sTree.substring(textPos, len(diresidues))
                            outputFile.write("C" + str(geneBoundaries[listIndex]['chromosome']) + ", " + outputTextPos + "\tPlus\t" + str(round(scores[child], 2)) + "\t" + sequence + "\t" + sequence + "\n")
                            
    urllib.urlopen("https://boglab.plp.iastate.edu/talent/jobcomplete/" + str(nodeID) + "/0")
