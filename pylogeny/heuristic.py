''' Define the interface for a heuristic in order to implement any manner of heuristic
for a combinatorial problem that can be abstracted into a state graph. In this case, a
phylogenetic tree space. '''

# Date:   Apr 1 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import os
from shutil import rmtree as removeFolder
from scoring import getLogLikelihood as ll
from executable import raxml
from tree import tree

# Constants

SMGR_MIN_NUM_LIKELIHOOD = 32

# Base Heuristic Superclass/Interface

class heuristic(object):

    ''' A base interface for a heuristic that explores a state graph. '''
    
    def __init__(self,G=None,start=None):
        self.landscape = G      # State Graph
        self.start     = start  # Start State
    def explore(self): pass
    def getStateGraph(self): return self.landscape
    def getStartState(self): return self.start

class phylogeneticLinearHeuristic(heuristic):
    
    ''' A base class for a heuristic that works on a phylogenetic landscape
    and only possesses a single path (of search). '''
    
    path, bestTree = list(), None
    def __init__(self,ls,startNode):
        super(phylogeneticLinearHeuristic,self
              ).__init__(ls,startNode)
    def getPath(self):     return self.path
    def getBestTree(self): return self.bestTree
    
# Heuristics

class parsimonyGreedy(phylogeneticLinearHeuristic):
    
    ''' Greedy (hill-climbing) landscape exploration by comparsion of parsimony. '''
    
    def __init__(self,ls,startNode):
        super(parsimonyGreedy,self).__init__(ls,startNode)
        
    def explore(self):
    
        ''' Perform greedy search of the landscape using
        a method of greed via parsimonious criterion. '''
        
        # Get starting tree, landscape.
        landscape = self.landscape
        ali       = landscape.alignment
        starting  = self.start
        cursor    = starting # Current tree.
        
        while (True):

            # Add this tree to the path.
            self.path.append(cursor)

            # Explore the current tree and only score via parsimony.
            landscape.exploreTree(cursor['index'])
            
            # Rank by parsimony.
            nodes = landscape.graph.neighbors(cursor['index'])
            nodes = [landscape.getNode(x) for x in nodes]
            nodes = [it for it in nodes if not it in self.path]
            nodes = sorted(nodes,key=lambda d: d['tree'].score[1])
            
            # Get best parsimony tree.
            best  = nodes[0]
            
            # Get parsimonies.
            bPars = best['tree'].score[1]
            cPars = cursor['tree'].score[1]
                 
            # Set new current tree if greater in parsimony.
            if (cPars > bPars): cursor = best
            else:
                self.bestTree = cursor
                break

class likelihoodGreedy(phylogeneticLinearHeuristic):
    
    ''' Greedy (hill-climbing) landscape exploration by comparsion of likelihood. '''
    
    def __init__(self,ls,startNode):
        super(likelihoodGreedy,self).__init__(ls,startNode)
        
    def explore(self):
    
        ''' Perform greedy search of the landscape using
        a method of greed via likelihood. '''
        
        # Get starting tree, landscape.
        landscape = self.landscape
        ali       = landscape.alignment
        starting  = self.start
        cursor    = starting # Current tree.
        
        while (True):

            # Add this tree to the path.
            self.path.append(cursor)

            # Explore the current tree and only score via parsimony.
            landscape.exploreTree(cursor['index'])
            
            # Rank by likelihood.
            nodes = landscape.getNeighborsFor(cursor['index'])
            nodes = [landscape.getNode(x) for x in nodes]
            nodes = [it for it in nodes if not it in self.path]
            map(lambda d: landscape.getVertex(d['index']).scoreLikelihood(),nodes)
            nodes = sorted(nodes,key=lambda d: d['tree'].score[0])
            
            # Get best likelihood tree.
            best  = nodes[-1]
            
            # Get likelihoods.
            bLikl = best['tree'].score[0]
            cLikl = cursor['tree'].score[0]
                 
            # Set new current tree if greater in likelihood.
            if (cLikl > bLikl): cursor = best
            else:
                self.bestTree = cursor
                break

class smoothGreedy(phylogeneticLinearHeuristic):
    
    ''' Parsimony-driven greedy landscape exploration 
    by comparsion of likelihoods. '''
    
    def __init__(self,ls,startNode):
        super(smoothGreedy,self).__init__(ls,startNode)
        
    def explore(self):
    
        ''' Perform greedy search of the landscape using
        a method of greed via parsimonious criterion and then
        performing final smoothing via likelihood on top 10% of
        1-SPR neighbors ranked on basis of parsimony. '''
        
        # Get starting tree, landscape.
        landscape = self.landscape
        ali       = landscape.alignment
        starting  = self.start
        cursor    = starting # Current tree.
        
        while (True):

            # Add this tree to the path.
            self.path.append(cursor)

            # Explore the current tree and only score via parsimony.
            landscape.exploreTree(cursor['index'])
            
            # Rank by parsimony.
            nodes = landscape.graph.neighbors(cursor['index'])
            nodes = [landscape.getNode(x) for x in nodes]
            nodes = [it for it in nodes if not it in self.path]
            nodes = sorted(nodes,key=lambda d: d['tree'].score[1])
            
            # Get best parsimony tree.
            best  = nodes[0]
            
            # Get parsimonies.
            bPars = best['tree'].score[1]
            cPars = cursor['tree'].score[1]
            
            # Get likelihood of current tree if necessary.
            cL = cursor['tree'].score[0]
            if (cL == None):
                cL = ll(cursor['tree'],ali)
                cursor['tree'].score = (cL,cPars)            
            
            # Set new current tree if greater in parsimony.
            if (cPars > bPars): cursor = best
            else: # Smoothing via likelihood.
                # Check likelihood of some number of neighbors.
                numtodo = max(SMGR_MIN_NUM_LIKELIHOOD,len(nodes)/10)
                numtodo = min(numtodo,len(nodes))
                for node in nodes[:numtodo]:                
                    # Get log-likelihood of this neighbor.
                    if (node['tree'].score[0] == None):
                        ml = ll(node['tree'],ali)
                        node['tree'].score = (ml,node['tree'].score[1])
                    else: ml = node['tree'].score[0]
                    if (cL < ml):
                        # Short circuit.
                        cursor = node
                        continue
                # Rank by likelihood.
                visited = sorted(self.path,key=lambda d: d['tree'].score[0])
                nodes   = sorted(nodes,key=lambda d: d['tree'].score[0])
                # Get the best likelihood tree.
                vbest = visited[-1]
                best  = nodes[-1]
                if (vbest['tree'].score[0] > best['tree'].score[0]):
                    self.bestTree = vbest
                    break
                else: cursor = best
                
class RAxMLIdentify(phylogeneticLinearHeuristic):

    ''' RAxML-driven landscape evaluation
    of intermediate checkpoint trees output from the 
    RAxML executable. '''
    
    def __init__(self,ls,startNode,workdir='.rxml'):
        super(RAxMLIdentify,self).__init__(ls,startNode)
        self.executable = None
        self.workdir    = workdir        
        self.startFile  = os.path.join(workdir,'start')
        self.logFile    = os.path.join(workdir,'logfi')
        self.__setupWorkDir__()
        self.__setupExecutable__()
    
    def __setupWorkDir__(self):
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)
        o = open(self.startFile,'w')
        o.write(self.start['tree'].newick)
        o.close()
    
    def __setupExecutable__(self):
        ali = self.landscape.alignment
        phy = ali.getPhylip()
        r   = raxml(phy,'rx',interTrees=True,startingTree=self.startFile,
                    log=self.logFile,wdir=self.workdir)
        self.executable = r
    
    def __readLogFile__(self):
        checklogf = os.path.join(self.workdir,'RAxML_log.rx')
        if not os.path.isfile(checklogf):
            raise IOError('RAxML output logfile not found (%s).' % (checklogf))
        o = open(checklogf,'r')
        lines = o.readlines()
        iters = []
        for it in lines: iters.append(it.strip('\n').split(' '))
        o.close()
        return iters
        
    def __readIterTrees__(self,iters):
        trees = []
        for it in iters:
            chfi       = 'RAxML_checkpoint.rx.%s' % (it[-1])
            checkpoint = os.path.join(self.workdir,chfi)
            o          = open(checkpoint)
            tr         = tree(o.read().strip('\n'))
            o.close()
            trees.append(tr)
        return trees
    
    def explore(self):
        
        # Get the landscape.
        landscape = self.landscape
        
        try:
            
            # Run the RAxML executable, read output.
            self.executable.run()
            checkpoints = self.__readLogFile__()
            
            # Get all intermediate trees.
            trees = self.__readIterTrees__(checkpoints)
            
            # See if these trees are already in the landscape.
            for t in trees:
                find = landscape.findTreeTopology(t.newick)
                if not find:
                    t.origin = 'RAxML'
                    find = landscape.addTree(t)
                self.path.append(landscape.getNode(find))
                if t == trees[-1]: self.bestTree = find
                
            # Remove work directory.
            removeFolder(self.workdir)
            
        except Exception,e:
            removeFolder(self.workdir)
            raise IOError('RAxML I/O failure for reason [%s].' % (str(e)))
            
        
