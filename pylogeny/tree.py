''' Container definition for (phylogenetic) bifurcating or multifurcating trees defined using Newick strings. '''

# Date:   Oct 20 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

# Imports

import newick, rearrangement
from math import factorial as fact

# Function Defitions

numberRootedTrees   = lambda t: numberUnrootedTrees(t+1)
numberUnrootedTrees = lambda t: (fact(2*(t-1)-3))/((2**(t-3))*fact(t-3))

# Class Definitions

class tree(object): # TODO: Integrate with P4 Tree class (?).
    
    ''' Defines a single (phylogenetic) tree by newick string;
    can possess other metadata. '''
    
    def __init__(self,newi='',check=False):
        
        ''' If enabled, "check" will force the structure to reroot
        the given Newick string tree to a lowest-order leaf in order
        to ensure a consistent Newick string among any duplicate
        topologies. '''
        
        self.name   = ''
        self.score  = None
        self.origin = None
        self.newick = newi
        if (check): self._checkNewick(newi)        
        else:       self.setNewick(newi)
    
    # Getters, Mutators
    
    def getName(self):     return self.name
    def setName(self,n):   self.name = n
    def getScore(self):    return self.score
    def setScore(self,s):  self.score = s
    
    def getOrigin(self):   return self.origin
    def setOrigin(self,o):
        
        ''' Set the "origin" or specification of where this tree
        was acquired or constructed from; a string. '''
        
        self.origin = o    
        
    def getNewick(self):   return self.newick
    def setNewick(self,n):
        
        ''' Set Newick string to n; also reacquires corresponding
        "structure" or Newick string without branch lengths. '''
        
        self.newick = n
        self.struct = self._getStructure()   
        
    def getStructure(self):
        
        ''' Returns "structure", a Newick string without branch lengths. '''
        return self.struct 
    
    def getSimpleNewick(self):

        ''' Return a Newick string with all taxa name replaced with
        successive integers. '''
    
        o = newick.parser(self.newick).parse()
        n = getAllLeaves(o)
        for _ in xrange(1,len(n)+1): n[_].label = str(_)
        return str(o) + ';'        
    
    # Other Functionality
    
    def toTopology(self):
        
        ''' Return a rearrangement.topology instance for this tree to allow for 
        rearrangement of the actual structure of the tree. '''
        
        t = rearrangement.topology()
        t.fromNewick(self.newick)
        return t
    
    # Internals
        
    def __eq__(self,o): return (o.struct == self.struct)
    def __ne__(self,o): return not (self.__eq__(o))
    def __str__(self): return self.newick
    
    def _checkNewick(self,newi):
        
        ''' PRIVATE: Run the Newick string through a 
        parse pass and reroot to lowest-order leaf in
        order to ensure a consistent Newick string. '''
        
        newi        = newi.strip('\n').strip(';') + ';'   
        prsd        = newick.parser(newi).parse()
        self.newick = rearrangement.topology(prsd).toNewick()
        self.struct = self._getStructure(prsd)
    
    def _getStructure(self,prsd=None):
        
        ''' PRIVATE: Acquires a newick string without any
        defined branch lengths. '''
        
        if prsd: p = prsd
        else: p = newick.parser(self.newick).parse()    
        removeBranchLengths(p)
        return str(p) + ';'  
    
class treeSet(object):
    
    ''' Represents an ordered, unorganized collection of trees 
    that do not necessarily comprise a combinatorial space. '''
    
    def __init__(self):
       
        self.trees = list()    
    
    def addTree(self,tr): 

        ''' Add a tree object to the collection. '''
    
        self.trees.append(tr)

    def addTreeByNewick(self,newick):
        
        ''' Add a tree to the structure by Newick string. '''
        
        t = tree(newick)
        self.addTree(t)    
        
    def removeTree(self,tr):
        
        ''' Remove a tree object from the collection if present. '''
        
        if tr in self.trees:
            self.trees.remove(tr)
    
    def indexOf(self,tr):
        
        ''' Acquire the index in this collection of a tree object. 
        Returns -1 if not found. '''
        
        if tr in self.trees: return self.trees.index(tr)
        else: return -1
    
    def __getitem__(self,i): return self.trees[i]
    
    def toTreeFile(self,fout):
        
        ''' Output this landscape as a series of trees, separated by
        newlines, as a text file saved at the given path. '''        
        
        o = open(fout,'w')
        o.write(str(self))
        o.close()
        return fout
        
    def __str__(self):
        return '\n'.join([t.getNewick() for t in self.trees])