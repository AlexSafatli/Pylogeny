''' Using the term borrowed from nomenclature of a bipartite graph, a bipartition for a phylogenetic tree coincides with the definition of two disjoint sets U and V . A branch in a phylogenetic tree defines a single bipartition that divides the tree into two disjoint sets U and V . The set U comprises all of the children leaf of the subtree associated with that branch. The set V contains the rest of the leaves or taxa in the tree. This package handle operations involving the representation of tree bipartitions. '''

# Date:   Apr 15 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import newick
from numpy import median

# Main Bipartition Object

class bipartition(object):

    ''' A tree bipartition. Requires a tree topology. '''

    def __init__(self,topol,bra=None):
        self.topology = topol
        self.branch   = bra
        self.btuple   = None # Tuple list representation.
        self.strrep   = None # String representation
        self.shortstr = ''   # Shorter string representation
        self.shortmap = None # Shorter string mapping of taxa to symbols
        self.reconfis = None # Possible reconfigurations as rearrangement objects.
        self._getStringRepresentation()
        self._getShortStringRepresentation()
        self._getBranchListRepresentation()
    
    def __hash__(self): return hash(self.shortstr)
    
    def __eq__(self,o):
        
        if (type(o) != bipartition): return False
        l,r   = self.strrep
        lo,ro = o.strrep
        return (lo==l and ro==r) or (ro==l and lo==r)
    
    def __ne__(self,o): return not self.__eq__(o)
    
    def _getStringRepresentation(self):
        
        ''' PRIVATE: Obtain string representation of a branch in a topology. '''
        
        if (self.branch == None): return
        if (self.topology == None):
            raise ValueError('Topology does not exist or is equal to None.')
        z            = self.topology.getStrBipartitionFromBranch(self.branch)
        l,r          = sorted(z[0]),sorted(z[1])
        self.strrep  = (l,r)
        
    def _getShortStringRepresentation(self):
        
        ''' PRIVATE: Obtain short string representation of branch/bipartition. '''
        
        ind  = 65 # 'A'
        l,r  = self.strrep
        x    = max((l,r),key=lambda d:len(d))
        if x == l: y = r
        else: y = l
        l,r  = x,y
        o    = ''
        lbls = [x for x in l] + [x for x in r]
        lblm = {}
        lbls = sorted(lbls) # Sort them.
        for i in lbls:
            lblm[i] = chr(ind)
            ind += 1
        for i in l: o += lblm[i]
        o += ':'
        for i in r: o += lblm[i]
        self.shortstr = o
        self.shortmap = lblm
        
    def _getBranchListRepresentation(self):
        
        l = [self.branch] + newick.getAllBranches(self.branch)
        r = [branch for branch in self.topology.getBranches(
            ) if not branch in l]
        self.btuple = (l,r)       
        
    def _getBranchFromString(self):
        
        ''' PRIVATE: Obtain branch information from a string representation. '''
        
        if not self.strrep or type(self.strrep) != tuple:
            raise IOError('Could not read string representation of bipartition.')
        self.branch = self.topology.getBranchFromStrBipartition(self.strrep)
        
    def _makeSPRRearrangements(self):
        
        ''' PRIVATE: Perform possible rearrangements. '''
        
        self.reconfis = self.topology.allSPRForBranch(self.branch)

    def fromStringRepresentation(self,st):
        
        ''' Acquire all component elements from a string representation 
        of a bipartition. '''
        
        self.strrep = st
        self._getBranchFromString()
        self._getShortStringRepresentation()
        
    def getBranch(self):
        
        ''' Get branch corresponding to this bipartition. '''
        
        return self.branch
    
    def getStringRepresentation(self):
        
        ''' Get the string representation corresponding to this bipartition. '''
        
        return self.strrep
    
    def getShortStringRepresentation(self):
        
        ''' Get the shorter string representation corresponding to this bipartition. '''
        
        return self.shortstr
    
    def getShortStringMappings(self):
        
        ''' Get the mapping of symbols from taxa names for the shorter string representation. '''
        
        return self.shortmap
    
    def getBranchListRepresentation(self):
        
        ''' Get the tuple of lists of branches that represent this bipartition. '''
        
        return self.btuple
    
    def getSPRRearrangements(self):
        
        ''' Return the set of all scores related to this bipartition. '''
        
        if (self.reconfis == None): self._makeSPRRearrangements()
        return self.reconfis
    
    def getSPRScores(self,ls,node=None):
        
        ''' Given a landscape, return all possible scores, not actively performing
        scoring if not done. '''
        
        # Get starting information.
        scores = list() # Output scoreset.
        if not node:
            origin = self.topology # Starting topology.
            if not hasattr(origin,'orig'): ornewi = origin.toNewick()
            else: ornewi = origin.orig
            # Acquire node corresponding to this topology.
            start = ls.findTreeTopology(ornewi)
            if (start == None): return list() # Could not find.
        else: start = node 
            
        # Get neighbors.
        neighbors = ls.graph.neighbors(start)
        
        # Investigate all rearrangements.
        rearrangements = self.getSPRRearrangements()
        resultants = [x.toTree().struct for x in rearrangements]
        for neighbor in neighbors:
            node   = ls.getNode(neighbor)
            tree   = ls.getTree(neighbor)
            struct = tree.struct
            if struct in resultants:
                scores.append(tree.score[0])
        return scores
    
    def getMedianSPRScore(self,ls,node=None):
        
        ''' Given a landscape, return the median SPR score. '''
        
        li = [x for x in self.getSPRScores(ls,node) if x]
        if len(li) > 0: return median(li)
        else: return None
        
    def getBestSPRScore(self,ls,node=None):
        
        ''' Given a landscape, return the best SPR score. '''
        
        li = [x for x in self.getSPRScores(ls,node) if x]
        if len(li) > 0: return max(li)
        else: return None