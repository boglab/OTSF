from cython.operator cimport dereference
from cython.operator import dereference

import numpy as np
cimport numpy as np

from libcpp.stack cimport stack
from libc.stdlib cimport malloc, free, calloc
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.map cimport map

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

# X replaces * in RVD names

#runtime scoring matrix type
scoringMatrixDtype = np.dtype([
    ('HD', 'float32'),
    ('NI', 'float32'),
    ('NG', 'float32'),
    ('NN', 'float32'),
    ('NS', 'float32'),
    ('NX', 'float32'),
    ('HG', 'float32'),
    ('HA', 'float32'),
    ('ND', 'float32'),
    ('NK', 'float32'),
    ('HI', 'float32'),
    ('HN', 'float32'),
    ('NA', 'float32'),
    ('IG', 'float32'),
    ('HX', 'float32'),
    ('SX', 'float32'),
    ('NH', 'float32'),
    ('YG', 'float32'),
    ('SN', 'float32'),
    ('SS', 'float32'),
    ('NC', 'float32'),
    ('HH', 'float32'),
    ('XX', 'float32'),
])

#compile time scoring matrix type
cdef packed struct scoringMatrixRecord:
    np.float32_t HD
    np.float32_t NI
    np.float32_t NG
    np.float32_t NN
    np.float32_t NS
    np.float32_t NX
    np.float32_t HG
    np.float32_t HA
    np.float32_t ND
    np.float32_t NK
    np.float32_t HI
    np.float32_t HN
    np.float32_t NA
    np.float32_t IG
    np.float32_t HX
    np.float32_t SX
    np.float32_t NH
    np.float32_t YG
    np.float32_t SN
    np.float32_t SS
    np.float32_t NC
    np.float32_t HH
    #other
    np.float32_t XX
        
ctypedef packed struct talentQueueItem:
    ulong nid
    double score

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
        
        
    cpdef ulong s(self, ulong v):
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




def ScoreTalentTask(querySequence, outputFilepath, bool revComp, geneBoundaries, np.ndarray[scoringMatrixRecord] scoringMatrix, PySSTree psTree):
    
    cdef SSTree *sTree = psTree.thisptr
    cdef double bestScore = 0
    
    #querySequence = querySequence.replace("*", "X").upper().rstrip()
    #diresidues = querySequence.split(' ')
    
    querySequence = querySequence.upper().rstrip()
    diresidues = querySequence.replace("*", "X").split()
    
    for diresidue in diresidues:
        bestScore += np.amin(scoringMatrix[diresidue])
    
    cdef double cutoffScore = 2.5 * bestScore
    
    with open(outputFilepath + ".txt", "w") as tabOutFile:
        with open(outputFilepath + ".gff3", "w") as gffOutFile:
            
            if not revComp:
                tabOutFile.write("table_ignores:Plus strand sequence\n")
            
            tabOutFile.write("Best Possible Score:" + str(round(bestScore, 2)) + "\n")
            tabOutFile.write('Genome Coordinates\tStrand\tScore\tTarget Sequence\tPlus strand sequence\n')
            
            gffOutFile.write("##gff-version 3\n")
            
            _ScoreTalentTask(querySequence, diresidues, len(diresidues), tabOutFile, gffOutFile, False, cutoffScore, baseMap, geneBoundaries, scoringMatrix, sTree)
            
            if revComp:
                _ScoreTalentTask(querySequence, diresidues, len(diresidues), tabOutFile, gffOutFile, True, cutoffScore, revBaseMap, geneBoundaries, scoringMatrix, sTree)
        

cdef char* reverseComplement(char* sequence, unsigned int sequenceLength):
    cdef char *new_sequence = <char*> calloc(sequenceLength, sizeof(char))
    cdef char base
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
    
cdef _ScoreTalentTask(querySequence, diresidues, unsigned int diresiduesLength, tabOutFile, gffOutFile, bool revComp, double cutoffScore, map[uchar, int] baseMap, geneBoundaries, np.ndarray[scoringMatrixRecord] scoringMatrix, SSTree *sTree):
    
    cdef ulong k, childNid, textPos
    cdef unsigned int parentDepth, depth
    cdef stack[talentQueueItem*] openSet
    cdef char startChar = 'T'
    cdef char *sequence, *revcomp_sequence 
    
    cdef double childScore
    cdef talentQueueItem* node
    cdef talentQueueItem* child
    
    cdef talentQueueItem *startNodeItem
    cdef ulong startNode
    cdef uchar edgeChar
    
    if revComp:
        
        diresidues.reverse()
        startNodeItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
        startNodeItem.nid = 0
        startNodeItem.score = 0
        
    
    else:
        startNodeItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
        startNodeItem.nid = sTree.search(<uchar*> &startChar, 1)
        startNodeItem.score = 0
        
        #startChar = 'A'
        #startNodeItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
        #startNodeItem.nid = sTree.search(<uchar*> &startChar, 1)
        #startNodeItem.score = 0
        #
        #startChar = 'C'
        #startNodeItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
        #startNodeItem.nid = sTree.search(<uchar*> &startChar, 1)
        #startNodeItem.score = 0
        #
        #startChar = 'G'
        #startNodeItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
        #startNodeItem.nid = sTree.search(<uchar*> &startChar, 1)
        #startNodeItem.score = 0
        
    
    openSet.push(startNodeItem)
    
    while not openSet.empty():
        
        node = openSet.top()
        openSet.pop()
        
        if not revComp:
            parentDepth = sTree.depth(node.nid) - 1
        else:
            parentDepth = sTree.depth(node.nid)
        
        if parentDepth < diresiduesLength:
        
            childNid = sTree.firstChild(node.nid)
            
            while childNid != 0:
                
                childScore = node.score
                
                k = 1
                depth = parentDepth
                edgeChar = sTree.edge(childNid, k)
                
                badChar = False
                
                while depth < diresiduesLength and edgeChar != '\x00':
                    
                    if edgeChar != '\x41' and edgeChar != '\x43' and edgeChar != '\x47' and edgeChar != '\x54':
                        badChar = True
                        break
                    
                    diresidue = diresidues[depth]
                    baseIndex = baseMap[edgeChar]
                    
                    if diresidue in scoringMatrixDtype.fields:
                        childScore += scoringMatrix[diresidue][baseIndex]
                    else:
                        childScore += scoringMatrix['XX'][baseIndex]
                    
                    k += 1
                    depth += 1
                    edgeChar = sTree.edge(childNid, k)
                
                if not badChar and childScore < cutoffScore:
                    
                    if revComp and depth == diresiduesLength:
                        
                        if edgeChar == '\x00' or edgeChar == '\x41':
                            
                            if edgeChar == '\x00':
                                childNid = sTree.child(childNid, 'A')
                            
                            
                            child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                            child.nid = childNid
                            child.score = childScore
                            openSet.push(child)
                            
                    else:
                        child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                        child.nid = childNid
                        child.score = childScore
                        openSet.push(child)
                        
                
                childNid = sTree.sibling(childNid)
                
        else:
            
            if not sTree.isleaf(node.nid):
                
                childNid = sTree.firstChild(node.nid)
            
                while childNid != 0:
                    
                    child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                    child.nid = childNid
                    child.score = node.score
                    openSet.push(child)
                    
                    childNid = sTree.sibling(childNid)
            
            else:
                
                if revComp:
                    textPos = sTree.textpos(node.nid)
                else:
                    textPos = sTree.textpos(node.nid) + 1
                    
                listIndex = geneBoundaries.bisect_right({'pos': textPos})
                
                if listIndex != len(geneBoundaries):
                    
                    sequence = <char*> sTree.substring(textPos, diresiduesLength)
                    
                    if listIndex == 0:
                        outputTextPos = textPos
                    else:
                        outputTextPos = textPos - geneBoundaries[listIndex - 1]['pos']
                    
                    
                    cname = str(geneBoundaries[listIndex]['chromosome'])
                    
                    outputScore = round(node.score, 2)
                    
                    if not revComp:
                        #tabOutFile.write("C" + cname + ", " + str(outputTextPos) + "\tPlus\t" + str(outputScore) + "\t" + sequence + "\t" + sequence + "\n")
                        tabOutFile.write("C%s, %lu\t%s\t%.2lf\t%s\t%s\n" % (cname, outputTextPos, "Plus", outputScore, sequence, sequence))
                        gffOutFile.write("%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;\n" % (cname, "TALESF", "TAL_effector_binding_site", outputTextPos, outputTextPos + diresiduesLength - 1, node.score, "+", querySequence, sequence))
                    else:
                        revcomp_sequence = reverseComplement(sequence, diresiduesLength)
                        tabOutFile.write("C%s, %lu\t%s\t%.2lf\t%s\t%s\n" % (cname, outputTextPos, "Minus", outputScore, revcomp_sequence, sequence))
                        gffOutFile.write("%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;\n" % (cname, "TALESF", "TAL_effector_binding_site", outputTextPos, outputTextPos + diresiduesLength - 1, node.score, "-", querySequence, revcomp_sequence))
                        free(revcomp_sequence)
                        
                    free(sequence)
                    
        free(node)
