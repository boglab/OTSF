from cython.operator cimport dereference
from cython.operator import dereference

import numpy as np
cimport numpy as np

from libcpp.deque cimport deque
from libc.stdlib cimport malloc, free
from libcpp cimport bool
from libcpp.pair cimport pair

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
        
def ScoreTalentTask(querySequence, outputFilepath, subMatrix, baseMap, geneBoundaries, np.ndarray[scoringMatrixRecord] scoringMatrix, PySSTree psTree):
    
    cdef ulong startNode, k, childNid, textPos
    cdef unsigned int parentDepth, depth
    cdef SSTree *sTree = psTree.thisptr
    cdef deque[talentQueueItem*] *openSet = new deque[talentQueueItem*]()
    cdef char* startLetter = "T"
    cdef double bestScore = 0
    cdef double childScore
    cdef talentQueueItem* node
    cdef talentQueueItem* child
    
    querySequence = querySequence.replace("*", "Z").upper().rstrip()
    diresidues = querySequence.split(' ')
    
    cdef unsigned int diresiduesLength = len(diresidues)
    
    
    for diresidue in diresidues:
        bestScore += np.amin(scoringMatrix[diresidue])
    
    cdef double cutoffScore = 2.5 * bestScore
    
    
    startNode = sTree.search(<uchar*> startLetter, 1)
    
    cdef talentQueueItem *testItem = <talentQueueItem*> malloc(sizeof(talentQueueItem))
    testItem.nid = startNode
    testItem.score = 0
    
    openSet.push_back(testItem)
    
    with open(outputFilepath, "w") as outputFile:
        
        #tableIgnores = []
        #if not revComp:
        
        #tableIgnores.append("Plus strand sequence")
        #if len(tableIgnores) > 0:
        outputFile.write("table_ignores:Plus strand sequence\n") #+ string.join(tableIgnores, ",") + "\n")
            
        outputFile.write("Best Possible Score:" + str(round(bestScore, 2)) + "\n")
        outputFile.write('Genome Coordinates\tStrand\tScore\tTarget Sequence\tPlus strand sequence\n')
        
        while not openSet.empty():
            
            node = openSet.back()
            openSet.pop_back()
            parentDepth = sTree.depth(node.nid) - 1
            
            if parentDepth < diresiduesLength:
            
                childNid = sTree.firstChild(node.nid)
                
                while childNid != 0:
                    
                    childScore = node.score
                    
                    k = 1
                    depth = parentDepth
                    edgeChar = "%c" % sTree.edge(childNid, k)
                    
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
                            childScore += scoringMatrix['ZZ'][baseIndex]
                        
                        k += 1
                        depth += 1
                        edgeChar = "%c" % sTree.edge(childNid, k)
                    
                    if not badChar and childScore < cutoffScore:
                        child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                        child.nid = childNid
                        child.score = childScore
                        openSet.push_back(child)
                    
                    childNid = sTree.sibling(childNid)
                    
            else:
                
                if not sTree.isleaf(node.nid):
                    
                    childNid = sTree.firstChild(node.nid)
                
                    while childNid != 0:
                        
                        child = <talentQueueItem*> malloc(sizeof(talentQueueItem))
                        child.nid = childNid
                        child.score = node.score
                        openSet.push_back(child)
                        
                        childNid = sTree.sibling(childNid)
                
                else:
                    
                    textPos = sTree.textpos(node.nid) + 1
                    listIndex = geneBoundaries.bisect_right({'pos': textPos})
                    
                    if listIndex != len(geneBoundaries):
                        
                        if listIndex == 0:
                            outputTextPos = str(textPos)
                        else:
                            outputTextPos = str(textPos - geneBoundaries[listIndex - 1]['pos'])
                        
                        sequence = <char*> sTree.substring(textPos, diresiduesLength)
                        outputFile.write("C" + str(geneBoundaries[listIndex]['chromosome']) + ", " + outputTextPos + "\tPlus\t" + str(round(node.score, 2)) + "\t" + sequence + "\t" + sequence + "\n")
            
            free(node)