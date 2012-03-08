from cython.operator cimport dereference, preincrement
from cython.operator import dereference, preincrement

import numpy as np
cimport numpy as np

from libcpp.stack cimport stack
from libc.stdlib cimport malloc, free, calloc
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string

from libc.math cimport round as cround
import time
#https://groups.google.com/d/topic/cython-users/4h2_AYi_ncA/discussion
cdef extern from * namespace "SSTree":
    enum io_action:
        nop, load_from, save_to
    
cdef extern from "SSTree.h":
    
    ctypedef unsigned char uchar
    ctypedef unsigned long ulong
    
    cdef cppclass SSTree:
        SSTree(uchar *text, ulong n)
        SSTree(uchar *text, ulong n, bool deletetext, unsigned samplerate, io_action IOaction, char* filename)
        ulong root()
        bool isleaf(ulong v)
        ulong child(ulong v, uchar c)
        ulong child(ulong v, uchar* s, ulong& l)
        ulong firstChild(ulong v)
        ulong sibling(ulong v)
        ulong parent(ulong v)
        uchar edge(ulong v, ulong d)
        char* edge(ulong v)
        
        char* edgePrefix(ulong, ulong)
        
        char* pathlabel(ulong v)
        
        char* pathlabel(ulong v, ulong l)
        
        char* substring(ulong i, ulong k)
        ulong depth(ulong v)
        ulong nodeDepth(ulong v)
        ulong lca(ulong v, ulong w)
        ulong lce(ulong i, ulong j)
        ulong sl(ulong v)
        ulong inorder(ulong v)
        
        ulong nextinpreorder(ulong)
        
        ulong rightmost(ulong v)
        ulong leftmost(ulong v)
        
        pair[ulong, ulong] subtreeInterval(ulong v)
        
        ulong leftrank(ulong v)
        ulong numberofnodes(ulong v)
        ulong numberofleaves(ulong v)
        ulong textpos(ulong v)
        ulong isOpen(ulong v)
        
        ulong SA(ulong)
        
        
        ulong search(uchar* pattern, ulong l)
        ulong lcaParen(ulong v, ulong w)
        void PrintTree(ulong v, ulong d)
        
        ulong LF(uchar c, ulong& sp, ulong& ep)
        ulong inverseSuffixLink(ulong nodeId, uchar c)

ctypedef bool (*talentQueueItemPFunctionPtr)(talentQueueItem*,talentQueueItem*)

cdef extern from "<algorithm>" namespace "std":
    void sort(vector[talentQueueItem*].iterator first, vector[talentQueueItem*].iterator last, talentQueueItemPFunctionPtr)
    void reverse(vector[int].iterator first, vector[int].iterator last)
    
cdef bool talentQueueItemPCompare(talentQueueItem* i,talentQueueItem* j):
    return i.score < j.score

# X replaces * in RVD names

cdef map[string, int] numpyMap
cdef char* tStr = 'HD'
numpyMap[string(tStr)] = 0
tStr = 'NI'
numpyMap[string(tStr)] = 1
tStr = 'NG'
numpyMap[string(tStr)] = 2
tStr = 'NN'
numpyMap[string(tStr)] = 3
tStr = 'NS'
numpyMap[string(tStr)] = 4
tStr = 'NX'
numpyMap[string(tStr)] = 5
tStr = 'HG'
numpyMap[string(tStr)] = 6
tStr = 'HA'
numpyMap[string(tStr)] = 7
tStr = 'ND'
numpyMap[string(tStr)] = 8
tStr = 'NK'
numpyMap[string(tStr)] = 9
tStr = 'HI'
numpyMap[string(tStr)] = 10
tStr = 'HN'
numpyMap[string(tStr)] = 11
tStr = 'NA'
numpyMap[string(tStr)] = 12
tStr = 'IG'
numpyMap[string(tStr)] = 13
tStr = 'HX'
numpyMap[string(tStr)] = 14
tStr = 'SX'
numpyMap[string(tStr)] = 15
tStr = 'NH'
numpyMap[string(tStr)] = 16
tStr = 'YG'
numpyMap[string(tStr)] = 17
tStr = 'SN'
numpyMap[string(tStr)] = 18
tStr = 'SS'
numpyMap[string(tStr)] = 19
tStr = 'NC'
numpyMap[string(tStr)] = 20
tStr = 'HH'
numpyMap[string(tStr)] = 21
tStr = 'XX'
numpyMap[string(tStr)] = 22
        
cdef map[uchar, int] baseMap

baseMap['A'] = 0
baseMap['C'] = 1
baseMap['G'] = 2
baseMap['T'] = 3

cdef map[uchar, int] revBaseMap

revBaseMap['A'] = 3
revBaseMap['C'] = 2
revBaseMap['G'] = 1
revBaseMap['T'] = 0

ctypedef packed struct talentQueueItem:
    ulong nid
    double score
    bool revComp


cdef class PySSTree:

    cdef SSTree *thisptr
    
    def __cinit__(self, char* text = NULL, char* action = NULL, char* fileName = NULL):
    
        cdef io_action realAction = nop
    
        if action is not NULL:
            if str(action) == "save":
                realAction = save_to
            elif str(action) == "load":
                realAction = load_from
            self.thisptr = new SSTree(<uchar*>text, len(text) + 1, False, 0, realAction, fileName)
        else:
            self.thisptr = new SSTree(<uchar*>text, len(text) + 1)


    def __dealloc__(self):
        del self.thisptr
        
        
    cpdef ulong root(self):
        """
        Returns the position of the root node in the parentheses sequence.
        """
        return self.thisptr.root()
        
        
    cpdef bool isleaf(self, ulong v):
        """
        Check if node v is a leaf node.

        Method doesn't check for an open parenthesis at position v (because we can assume that v is a node).
        """
        return self.thisptr.isleaf(v)
        
        
    cpdef firstChild(self, ulong v):
        """
        Returns the first child of the (internal) node v.
        """  
        return self.thisptr.firstChild(v)
    
    
    def children(self, ulong v):
        """
        Generator for all children of the (internal) node v.
        """
        cdef ulong child = self.thisptr.firstChild(v)
        
        while child != 0:
            yield child
            child = self.thisptr.sibling(child)
        
        
    cpdef ulong sibling(self, ulong v):
        """
        Returns the next sibling of the node v.
        """  
        return self.thisptr.sibling(v)

        
    cpdef ulong parent(self, ulong v):
        """
        Returns the parent of the node v.
        """  
        return self.thisptr.parent(v)
        
        
    cpdef edge(self, ulong v):
        """
        Returns the edge label of the incoming edge of the node v.
        """  
        return <char*> self.thisptr.edge(v)
        
        
    cpdef pathlabel(self, ulong v):
        """
        Returns the path label from root to the node v.
        """  
        return <char*> self.thisptr.pathlabel(v)
        
        
    cpdef substring(self, ulong i, ulong k):
        """
        Returns the substring of the original text sequence with that starts at position i and has length k.
        """  
        return <char*> self.thisptr.substring(i, k)
        
        
    cpdef ulong depth(self, ulong v):
        """
        Returns the string depth of the node v.
        """  
        return self.thisptr.depth(v)
        
        
    cpdef ulong nodeDepth(self, ulong v):
        """
        Returns the node depth of the node v.
        
        Node depth of root is 0.
        """  
        return self.thisptr.nodeDepth(v)
        
        
    cpdef ulong sl(self, ulong v):
        """
        Suffix link for internal nodes
        """  
        return self.thisptr.sl(v)
        
    cpdef ulong wl(self, ulong v, char* c):
        """
        Weiner link for internal nodes
        """  
        return self.thisptr.inverseSuffixLink(<uchar> dereference(c), v)
        
    cpdef ulong inorder(self, ulong v):
        """
        Returns the Inorder number of a node
        """  
        return self.thisptr.inorder(v)
        
        
    cpdef ulong rightmost(self, ulong v):
        """
        Returns the Right most leaf of the (internal) node v.
        """  
        return self.thisptr.rightmost(v)
        
        
    cpdef ulong leftmost(self, ulong v):
        """
        Returns the Left most leaf of the (internal) node v.
        """  
        return self.thisptr.leftmost(v)
        
        
    cpdef ulong leftrank(self, ulong v):
        """
        Returns the Left rank of a (leaf) node.
        """  
        return self.thisptr.leftrank(v)
        
        
    cpdef ulong numberofnodes(self, ulong v):
        """
        Returns the Number of nodes in a subtree rooted at the node v.
        """  
        return self.thisptr.numberofnodes(v)
        
        
    cpdef ulong numberofleaves(self, ulong v):
        """
        Returns the Number of leaf nodes in a subtree rooted at the node v.
        """  
        return self.thisptr.numberofleaves(v)
        
        
    cpdef ulong textpos(self, ulong v):
        """
        Returns the Suffix array value of the leaf node v.
        """  
        return self.thisptr.textpos(v)
        
        
    cpdef bool isOpen(self, ulong v):
        """
        Check for an open parentheses at index v.
        """  
        return self.thisptr.isOpen(v)
        
    def edgeChar(self, ulong v, ulong d):
        """
        Returns the d-th character of the label of the node v's incoming edge.
        """  
        return "%c" % self.thisptr.edge(v, d)
        
    cpdef ulong lca(self, ulong v, ulong w):
        """
        Returns the Lowest common ancestor of nodes v and w.
        """  
        return self.thisptr.lca(v, w)
        
    cpdef ulong lce(self, ulong i, ulong j):
        """
        Returns the longest common extension of text positions i and j.
        """  
        return self.thisptr.lce(i, j)
        
    cpdef ulong lcaParen(self, ulong v, ulong w):
        """
        Returns the lowest common ancestor of nodes v and w.
        
        Linear time solution for debugging.
        """  
        return self.thisptr.lcaParen(v, w)
        
    cpdef ulong search(self, char* pattern):
        """
        Search for a given pattern of length l from the suffix array.
        
        Returns a node from the suffix tree that is LCA of the matching leaf interval.
        """  
        return self.thisptr.search(<uchar*>pattern, len(pattern))
        
    cpdef ulong child(self, ulong v, char* c):
        """
        Returns the child node of v so that edge leading to the child node starts with the symbol c.
        
        Returns 0 if no such child exists
        """  
        return self.thisptr.child(v, <uchar> dereference(c))
        
        
    #def child(self, ulong v, char* s, ulong& l ):
    #    """
    #    Returns the child node of v so that edge leading to the child node (or its prefix) is the string s.
    #
    #    The parameter s must be '\0' terminated. 
    #    The last parameter stores the number of matched letters.
    #    Note that end marker '\0' will not be matched.
    #    FIXME This is a slow version
    #    """
    #    
    #    return self.thisptr.child(v, <uchar*>s, l)
        
        
    cpdef PrintTree(self, ulong v, ulong d):
        self.thisptr.PrintTree(v, d)


def ScoreTalentTask(querySequence, outputFilepath, bool revComp, geneBoundaries, np.ndarray[np.float32_t, ndim=2] scoringMatrix, PySSTree psTree):
    
    cdef SSTree *sTree = psTree.thisptr
    cdef double bestScore = 0
    
    querySequence = querySequence.upper().rstrip()
    py_diresidues = querySequence.replace("*", "X").split()
    
    cdef vector[int] diresidues
    cdef char* diresidue_c_str
    cdef string diresidue_string
    cdef int numpyIndex

    cdef char* unknown_diresidue_c_str = "XX"
    cdef int unknown_numpy_index = numpyMap[string(unknown_diresidue_c_str)]
    
    for diresidue in py_diresidues:

        diresidue_c_str = diresidue
        diresidue_string = string(diresidue_c_str)
        
        if numpyMap.count(diresidue_string):
            numpyIndex = numpyMap[diresidue_string]
        else:
            numpyIndex = unknown_numpy_index
        
        diresidues.push_back(numpyIndex)
        bestScore += np.amin(scoringMatrix[numpyIndex])
    
    cdef double cutoffScore = 3.0 * bestScore
    
    cdef vector[talentQueueItem*] results
    
    _ScoreTalentTask(querySequence, &diresidues, &results, False, cutoffScore, baseMap, scoringMatrix, sTree)
    
    if revComp:
        _ScoreTalentTask(querySequence, &diresidues, &results, True, cutoffScore, revBaseMap, scoringMatrix, sTree)
        
    _PrintTaskResults(querySequence, diresidues.size(), outputFilepath, &results, revComp, bestScore, geneBoundaries, sTree)

cdef char* reverseComplement(char* sequence, unsigned int sequenceLength):
    cdef char *new_sequence = <char*> calloc(sequenceLength, sizeof(char))
    cdef char base
    cdef unsigned int i
    
    for i in range(sequenceLength):
        base = sequence[sequenceLength - i - 1]
        if base == 'A':
            new_sequence[i] = 'T'
        elif base == 'C':
            new_sequence[i] = 'G'
        elif base == 'G':
            new_sequence[i] = 'C'
        elif base == 'T':
            new_sequence[i] = 'A'
    
    return new_sequence
    
cdef _ScoreTalentTask(querySequence, vector[int] *diresidues, vector[talentQueueItem*] *results, bool revComp, double cutoffScore, map[uchar, int] baseMap, np.ndarray[np.float32_t, ndim=2] scoringMatrix, SSTree *sTree):
    
    cdef ulong k, childNid, revCompChild
    cdef unsigned int parentDepth, depth
    cdef stack[talentQueueItem*] openSet
    cdef char startChar = 'T'
    
    cdef double childScore
    cdef talentQueueItem *node
    cdef talentQueueItem *child
    
    cdef talentQueueItem *startNodeItem
    cdef ulong startNode
    cdef uchar edgeChar
    
    cdef bool nodeOutput
    
    cdef int baseIndex, numpyIndex, listIndex, diresidue
    
    startNodeItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
    startNodeItem.score = 0
    
    if revComp:
        
        diresidues = new vector[int](dereference(diresidues))
        reverse(diresidues.begin(), diresidues.end())
        startNodeItem.nid = 0
    
    else:
        
        startNodeItem.nid = sTree.search(<uchar*> &startChar, 1)
    
    openSet.push(startNodeItem)
    
    while not openSet.empty():
        
        node = openSet.top()
        openSet.pop()
        
        nodeOutput = False
        
        if not revComp:
            parentDepth = sTree.depth(node.nid) - 1
        else:
            parentDepth = sTree.depth(node.nid)
        
        if parentDepth < diresidues.size():
        
            childNid = sTree.firstChild(node.nid)
            
            while childNid != 0:
                
                childScore = node.score
                
                k = 1
                depth = parentDepth
                edgeChar = sTree.edge(childNid, k)
                
                while depth < diresidues.size() and childScore <= cutoffScore and edgeChar != '\x00':
                    
                    if edgeChar != '\x41' and edgeChar != '\x43' and edgeChar != '\x47' and edgeChar != '\x54':
                        childScore = cutoffScore + 1
                        break
                    
                    diresidue = diresidues.at(depth)
                    baseIndex = baseMap[edgeChar]
                    
                    childScore += scoringMatrix[diresidue, baseIndex]
                        
                    k += 1                    
                    depth += 1

                    edgeChar = sTree.edge(childNid, k)
                
                if childScore <= cutoffScore:
                    
                    if revComp and depth == diresidues.size() and edgeChar != '\x00':
                        
                        if edgeChar == '\x41':
                                
                            child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                            child.nid = childNid
                            child.score = childScore
                            child.revComp = revComp
                            openSet.push(child)
                        
                    else:
                        child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                        child.nid = childNid
                        child.score = childScore
                        child.revComp = revComp
                        openSet.push(child)
                
                childNid = sTree.sibling(childNid)
                
        elif revComp and parentDepth == diresidues.size():
            
            revCompChild = sTree.child(node.nid, 'A')
            
            if revCompChild != 0:
                child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                child.nid = revCompChild
                child.score = node.score
                child.revComp = revComp
                openSet.push(child)
                
        else:
            
            if not sTree.isleaf(node.nid):
                
                childNid = sTree.firstChild(node.nid)
            
                while childNid != 0:
                    
                    child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                        
                    child.nid = childNid
                    child.score = node.score
                    child.revComp = revComp
                    openSet.push(child)
                    childNid = sTree.sibling(childNid)
            
            else:
                
                nodeOutput = True
                results.push_back(node)
        
        if not nodeOutput:   
            free(node)
            
    if revComp:
        free(diresidues)

cdef _PrintTaskResults(querySequence, unsigned int diresiduesLength, outputFilepath, vector[talentQueueItem*] *results, bool revComp, double bestScore, geneBoundaries, SSTree *sTree):
    
    sort(results.begin(), results.end(), talentQueueItemPCompare)
    
    cdef vector[talentQueueItem*].iterator it = results.begin()
    cdef talentQueueItem* node
    cdef int counter = 1
    cdef char *sequence, *revcomp_sequence
    cdef ulong textPos
    
    with open(outputFilepath + ".txt", "w") as tabOutFile:
        with open(outputFilepath + ".gff3", "w") as gffOutFile:
            with open(outputFilepath + ".bed", "w") as bedOutFile:
            
                gffOutFile.write("##gff-version 3\n")
                
                if not revComp:
                    tabOutFile.write("table_ignores:Plus strand sequence\n")
                    gffOutFile.write("#table_display_tags:target_sequence\n")
                else:
                    gffOutFile.write("#table_display_tags:target_sequence,plus_strand_sequence\n")
                    
                outputQuerySequence = querySequence.replace(" ", "_")
                    
                bedOutFile.write('track name="TAL Targets" description="Targets for RVD sequence ' + querySequence + '" visibility=2 useScore=1\n')
                
                roundedBestScore = str(round(bestScore, 2))
                gffOutFile.write("#Best Possible Score:" + roundedBestScore + "\n")            
                
                tabOutFile.write("Best Possible Score:" + roundedBestScore + "\n")
                tabOutFile.write('Genome Coordinates\tStrand\tScore\tTarget Sequence\tPlus strand sequence\n')
                
                while it != results.end():
                    
                    node = dereference(it)
                    
                    if node.revComp:
                        textPos = sTree.textpos(node.nid)
                    else:
                        textPos = sTree.textpos(node.nid) + 1
                        
                    listIndex = geneBoundaries.bisect_right({'pos': textPos})
                    
                    if listIndex != len(geneBoundaries):
                        
                        sequence = <char*> sTree.substring(textPos, diresiduesLength)
                        
                        if listIndex == 0:
                            outputTextPos = textPos + 1
                        else:
                            outputTextPos = textPos - geneBoundaries[listIndex - 1]['pos']
                            
                        cname = str(geneBoundaries[listIndex]['chromosome'])
                        
                        outputScore = round(node.score, 2)
                        
                        if not node.revComp:
                            tabOutFile.write("C%s, %lu\t%s\t%.2lf\t%s\t%s\n" % (cname, outputTextPos, "Plus", outputScore, sequence, sequence))
                            gffOutFile.write("Chr%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;plus_strand_sequence=%s;\n" % (cname, "TALESF", "TAL_effector_binding_site", outputTextPos, outputTextPos + diresiduesLength - 1, node.score, "+", outputQuerySequence, sequence, sequence))
                            # 0 based
                            bedOutFile.write("%s\t%d\t%d\tsite%d\t%.2lf\t%s\n" %  (cname, outputTextPos - 1, outputTextPos + diresiduesLength - 2, counter, cround(bestScore / node.score * 1000), '+'))
                        else:
                            revcomp_sequence = reverseComplement(sequence, diresiduesLength)
                            tabOutFile.write("C%s, %lu\t%s\t%.2lf\t%s\t%s\n" % (cname, outputTextPos, "Minus", outputScore, revcomp_sequence, sequence))
                            gffOutFile.write("Chr%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;plus_strand_sequence=%s;\n" % (cname, "TALESF", "TAL_effector_binding_site", outputTextPos, outputTextPos + diresiduesLength - 1, node.score, "-", outputQuerySequence, revcomp_sequence, sequence))
                            # 0 based
                            bedOutFile.write("%s\t%d\t%d\tsite%d\t%.2lf\t%s\n" %  (cname, outputTextPos - 1, outputTextPos + diresiduesLength - 2, counter, cround(bestScore / node.score * 1000), '-'))
                        
                        if revcomp_sequence != NULL:
                            free(revcomp_sequence)
                            revcomp_sequence = NULL
                            
                        free(sequence)
                        
                        counter += 1
                    
                    free(node)
                    preincrement(it)