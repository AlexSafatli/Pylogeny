''' Encapsulate a phylogenetic tree space. A phylogenetic landscape or tree
space refers to the entire combinatorial space comprising all possible
phylogenetic tree topologies for a set of n taxa. The landscape of n taxa can be
defined as consisting of a finite set T of tree topologies. Tree topologies can
be associated with a fitness function f(t_i) describing their fit. This forms a
discrete solution search space and finite graph (T, E) = G. E(G) refers to the
neighborhood relation on T(G). Edges in this graph are bidirectional and
represent transformation from one tree topology to another by a tree
rearrangement operator. An edge between t_i and t_j would be notated as e_{ij}
in E(G). '''

# Date:   Jan 24 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

# Imports

import networkx
import tree
import alignment
import base
from random import choice
from scoring import getParsimonyFromProfiles as parsimony, getLogLikelihood as\
     ll
from parsimony import profile_set as profiles
from networkx import components as comp, algorithms as alg
from base import patriciaTree
from tree import treeSet, numberRootedTrees, numberUnrootedTrees
from newick import newickParser, removeBranchLengths
from rearrangement import TYPE_NNI, TYPE_SPR, TYPE_TBR
postOrderTraversal = base.treeStructure.postOrderTraversal

LS_NOT_DEFINED = -1

# Graph Object

class graph(object):
    
    ''' Define an empty graph object. '''
    
    def __init__(self,gr=None,defWeight=0.):
        
        ''' Instantiate a graph. Default edge weights are 0.
        
        :param gr: a networkx graph object, if already exists. 
        :type gr: a :class:`networkx.Graph` object
        :param defWeight: the default edge weight of weights
        :type defWeight: a floating point number
        
        '''
        
        if gr == None: self.graph = networkx.Graph()
        else: self.graph = gr
        self.defaultWeight = defWeight
    
    def __len__(self):
        
        return len(self.graph.node)
    
    def __iter__(self):
        
        for node in self.graph.node.keys(): yield node

    def getSize(self):

        ''' Return the number of nodes in the graph.
        
        :return: an integer
        
        '''
        
        return len(self.graph.node)

    def getNodeNames(self):
        
        ''' Return the names of nodes in the graph.
        
        :return: a list of strings
        
        '''
        
        return self.graph.node.keys()
    
    def iterNodes(self):
        
        ''' Iterate over all node keys.
        
        :return: a generator of node keys
        
        '''
        
        for node in self.graph.node: yield node
        
    def getNodes(self):
        
        ''' Get all node values. Recommended to use an iterator for large
        graphs (iterNodes()).
        
        :return: a list of node values (whatever is associated with nodes)
        
        '''
        
        return self.graph.node.values()
    
    def getEdges(self):
        
        ''' Get all edges (as defined for NetworkX graphs). Recommended to use
        an iterator for large graphs.
        
        :return: a list of edges (and associated data)
        
        '''
        
        return [self.getEdge(i,j) for i,j in self.graph.edges_iter()]
    
    def getEdgesFor(self,i):
        
        ''' Get all edges associated with a certain node. 
        
        :param i: a node name
        :type i: a string
        :return: a list of edges (and associated data) for all neighbors
        
        '''
        
        return [self.getEdge(i,j) for j in self.graph.neighbors(i)]
    
    def getNode(self,i):
        
        ''' Get a single node by name. 
        
        :param i: a node name:
        :type i: a string
        :return: a dictionary (with node information)
        
        '''
        
        return self.graph.node[i] 
    
    def getEdge(self,i,j):
        
        ''' Get the data associated with an edge (including weight).
        
        :param i: a node name
        :type i: a string
        :param j: a node name
        :type j: a string
        :return: an edge (and associated data)
        
        '''
        
        return self.graph.get_edge_data(i,j)
    
    def isEdge(self,i,j):
        
        ''' See if an edge exists between two nodes.
        
        :param i: a node name
        :type i: a string
        :param j: a node name
        :type j: a string
        :return: a boolean
        
        '''
        
        return (self.getEdge(i,j) != None)
    
    def getNeighborsFor(self,i):
        
        ''' Get a list of all node names neighbor to a node. 
        
        :param i: a node name
        :type i: a string
        :return: a list of strings (node names)
        
        '''
        
        return self.graph.neighbors(i)  
    
    def getDegreeFor(self,i):

        ''' Return in- and out-degree for node named i.
        
        :param i: a node name
        :type i: a string
        :return: an integer
        
        '''
        
        return len(self.getNeighborsFor(i))
        
    def setDefaultWeight(self,w):
        
        ''' Set the default weight of edges (weight of edges if not overridden).
        
        :param w: a weight
        :type w: a floating point
        
        '''
        
        self.defaultWeight = float(w)
        
    def clearEdgeWeights(self):
        
        ''' Set all edge weights to the default edge weight. '''
        
        for edge in self.getEdges():
            edge['weight'] = self.defaultWeight
        
    def getNumComponents(self):
        
        ''' Get the number of components of the graph.
        
        :return: an integer
        
        '''
        
        return comp.number_connected_components(self.graph)
    
    def getComponents(self):
        
        ''' Get the connected components in the graph. '''
        
        return comp.connected_components(self.graph)
        
    def getComponentOfNode(self,i):
        
        ''' Get the graph component of a given node. '''
        
        return comp.node_connected_component(self.graph,i)

    def getCliques(self):
        
        ''' Get the cliques present in the graph. '''
        
        return alg.clique.find_cliques(self.graph)
    
    def getCliqueNumber(self):

        ''' Get the clique number of the graph.
        
        :return: an integer
        
        '''
        
        return alg.clique.graph_clique_number(self.graph)
    
    def getNumCliques(self):
        
        ''' Get the number of cliques found in the graph.
        
        :return: an integer
        
        '''
        
        return alg.clique.number_of_cliques(self.graph)
    
    def getCliquesOfNode(self,i):
        
        ''' Get the clique that a node corresponds to. '''
        
        return alg.clique.cliques_containing_node(self.graph,i)
    
    def getCenter(self):
        
        ''' Get the centre of the graph. '''
        
        return alg.distance_measures.center(self.graph)
    
    def getDiameter(self):
        
        ''' Acquire the diameter of the graph. '''
        
        return alg.distance_measures.diameter(self.graph)    
    
    def getMST(self):

        ''' Acquire the minimum spanning tree for the graph. '''
        
        return alg.minimum_spanning_tree(self.graph)
    
    def hasPath(self,nodA,nodB):

        ''' See if a path exists between two nodes. '''
        
        return alg.has_path(self.graph,nodA,nodB)

    def getShortestPath(self,nodA,nodB):
        
        ''' Get the shortest path between two nodes. '''
        
        return alg.shortest_path(self.graph,nodA,nodB)

    def getShortestPathLength(self,nodA,nodB):
        
        ''' Get the shortest path length between two nodes. '''
        
        return alg.shortest_path_length(self.graph,nodA,nodB)


class landscape(graph,treeSet):
    
    ''' Defines an entire phylogenetic tree space. '''
    
    def __init__(self,ali,starting_tree=None,root=True,operator='SPR'):
        
        ''' Initialize the landscape.

        :param ali: an alignment
        :type ali: an :class:`.alignment.alignment` object
        :param starting_tree: an optional tree object to start with
        :type starting_tree: a :class:`.tree.tree` object
        :param root: whether or not to compute an approximate\
        maximum likelihood tree (FastTree) or start the landscape\
        with a given starting tree.
        :type root: a boolean
        :param operator: a string that describes what operator the\
        landscape is mostly comprised of.
        :type operator: a string
        
        '''        
        
        super(landscape,self).__init__()
        
        # Fields
        self.alignment          = ali
        self.locks              = list()
        self.nextTree           = 0
        self.leaves             = None
        self.root               = None
        self.operator           = operator
        self.parsimony_profiles = None
        self.newickSearchDict   = dict()
        
        # Analyze alignment.
        if ali:
            # Get number of leaves.
            self.leaves = ali.getNumSeqs()
            # Get parsimony profile.
            self.parsimony_profiles = profiles(ali)           
                
        # Set up the root.
        if root:
            if starting_tree == None:
                approx      = ali.getApproxMLTree()
                self.root   = self.addTree(approx)
            else: self.root = self.addTree(starting_tree)
            # Score this tree.
            node  = self.getNode(self.root)
            tre   = node['tree']
            topol = tre.toTopology()
            new   = topol.toNewick()
            if not (tre.score):
                sc = parsimony(new,p)
                tre.score = (None,sc)            
        
    def getAlignment(self): 
        
        ''' Acquire the alignment object associated with this space. 
        
        :return: an :class:`.alignment.alignment` object
        
        '''
        
        return self.alignment
    
    def getNumberTaxa(self):
        
        ''' Return the number of different taxa present in any respective tree
        in the landscape. 
        
        :return: an integer
        
        '''
        
        return self.leaves

    def getPossibleNumberRootedTrees(self):
        
        ''' Assuming all of the trees in the space are rooted, return the
        maximum possible number of unrooted trees that can possibly be generated
        for the number of taxa of trees in the landscape.
        
        :return: an integer
        
        '''
        
        return numberRootedTrees(self.leaves)
    
    def getPossibleNumberUnrootedTrees(self):
        
        ''' Assuming all of the trees in the space are unrooted, return the 
        maximum possible number of unrooted trees that can possibly be generated
        for the number of taxa of trees in the landscape.
        
        :return: an integer
        
        '''
        
        return numberUnrootedTrees(self.leaves)
    
    def getRoot(self):
        
        ''' Returns the index to the root (starting) tree of the space.
        
        :return: an integer
        
        '''
        
        return self.root
    
    def getRootNode(self):
        
        ''' Returns the root (starting) tree of the space in its node form. 
        
        :return: a dictionary (with node information)
        
        '''
    
    def getRootTree(self):
        
        ''' Acquire the first tree that was placed in this space.
        
        :return: a :class:`.tree.tree` object
        
        '''
        
        return self.getTree(self.root)
    
    def setAlignment(self,ali):
        
        ''' Set the alignment present in this landscape. WARNING; will not
        modify existing scores. 
        
        :param ali: an alignment
        :type ali: an :class:`.alignment.alignment` object
        
        '''
        
        self.alignment          = ali
        self.leaves             = ali.getNumSeqs()
        self.parsimony_profiles = profiles(ali)
        
    def setOperator(self,op):
        
        ''' Set the operator assigned to this landscape.
        
        :param op: an operator (string description)
        :type op: a string
        
        '''
        
        self.operator = op
    
    # Node Management
                
    def _newNode(self,tobj,score=False):
        
        ''' PRIVATE: Add a new node. '''
        
        # Extract data from the object.
        if type(tobj) != tree.tree:
            raise TypeError(
                'Nodes in space must be constructed from tree objects.')
        name      = tobj.getName()
        structure = tobj.getStructure()
        
        # Add its structure to a dictionary structure.
        query = self.findTreeTopologyByStructure(structure)
        if query != None:
            raise AssertionError('Tree (%s) <%s> already exists in space!' % (
                str(name),str(structure)))
        i = len(self) # Get next possible value for insertion unique integer.
        self.newickSearchDict[structure] = i
        
        # Create the node.
        self.graph.add_node(i)
        node = self.graph.node[i]
        node['index']    = i
        node['explored'] = False
        node['tree']     = tobj
        node['failed']   = False
        
        # Preliminary scoring.
        if score:
            if tobj.score == None:
                # Get the parsimony since this is fast; set likelihood
                # to nothing.
                if self.alignment:
                    tobj.score = (None,parsimony(
                        tobj.newick,self.parsimony_profiles))
                else: tobj.score = (None,None)
            elif tobj.score[1] == None:
                # Get the parsimony since this is fast; set likelihood
                # to nothing.
                if self.alignment:
                    tobj.score = (tobj.score[0],parsimony(
                        tobj.newick,self.parsimony_profiles))  
        
        # Return the index.
        return i
    
    def getTree(self,i):
        
        ''' Get the object for a tree by its name.
        
        :param i: a tree name (usually an integer)
        :return: a :class:`.tree.tree` object
        
        '''
        
        if not i in self.graph.node: return None
        return self.getNode(i)['tree']
    
    def __iter__(self):
        
        for t in self.graph.nodes():
            yield self.getTree(t)

    def getVertex(self,i):
        
        ''' Acquire a vertex object from the landscape; this is a 
        high-level representation of a tree in the landscape with
        additional functionality. Object created upon invocation of
        this function.
        
        :param i: a tree name (usually an integer)
        :return: a :class:`.vertex` object
        
        '''
        
        return vertex(self.getNode(i),self)

    def removeTreeByIndex(self,i):
        
        ''' Remove a tree from the landscape by index.
        
        :param i: a tree name (usually an integer)
        :return: a boolean (success or failure)
        
        '''
        
        tr = self.getTree(i)
        if (tr == None): return False
        self.graph.remove_node(i)
        del self.newickSearchDict[tr.getStructure()]
        return True

    def removeTree(self,tree):
        
        ''' Remove a tree from the landscape by object.
        
        :param tree: a tree that exists in the landscape
        :type tree: a :class:`.tree.tree` object
        :return: a boolean (success or failure)
        
        '''

        t = self.indexOf(tree)
        if (t < 0): return False
        return self.removeTreeByIndex(t)
        
    def addTreeByNewick(self,newick,score=True,check=True,struct=None):
        
        ''' Add tree to the landscape by Newick string. Will return index.
        
        :param newick: a Newick string
        :type newick: a string
        :param score: defaults to True, whether to score this tree or not
        :type score: a boolean
        :return: the index of the tree
        
        '''
        
        return self.addTree(None,score=score,check=check,newick=newick,
                            struct=struct)
        
    def addTree(self,tr,score=True,check=True,newick=None,struct=None):
        
        ''' Add a tree to the landscape. Will return its index. 

        :param tr: a tree
        :type tr: a :class:`.tree.tree` object
        :param score: defaults to True, whether to score this tree or not
        :type score: a boolean
        :return: the index of the tree
        
        '''
            
        # (Re)construct the tree ensuring consistent configuration/naming.
        if tr == None and newick == None:
            raise ValueError(
                'Require tree object OR Newick string for tree addition.')
        if newick == None: newick = tr.toNewick()
        tobj = tree.tree(newick,check=check,structure=struct)
        
        # See if needs to be scored.
        if not score and tr != None: tobj.setScore(tr.getScore())
        
        # Origin needs to be defined?
        if tr != None: tobj.setOrigin(tr.getOrigin())
        
        # Add to graph.
        index = self._newNode(tobj,score=score)
        if (self.root is None):
            if (len(self) != 1):
                raise AssertionError('Landscape has no root with 1+ tree.')
            self.root = index
        return index
    
    def exploreRandomTree(self,i,type=TYPE_SPR):
        
        ''' Acquire a single neighbor to a tree in the landscape by performing a
        random rearrangement of type SPR (by default), NNI, or TBR -- this is
        done by performing a rearrangement on a random branch in the topology.
        Rearrangement type is provided as a rearrangement module type definition
        of form, for example, TYPE_SPR, TYPE_NNI, etc.
        
        :param i: a tree index
        :param type: the type of rearrangement (e.g., TYPE_SPR, TYPE_NNI)
        :return: the new tree index or None in case of failure
        
        '''
        
        # Get parsimony profiles.
        p = self.parsimony_profiles
        
        # Check node.
        node  = self.getNode(i)
        if (node['explored']): return None
        tre   = node['tree']
        topol = tre.toTopology()
        new   = topol.toNewick()
        isvi  = self._isViolating(topol)
        
        # Check for locks.
        if hasattr(self,'locks'):
            for lock in self.locks:
                hasthis = topol.getBranchFromBipartition(lock)
                if (hasthis): topol.lockBranch(hasthis)
        
        # Perform an exploration (not yet done).
        branches = topol.getBranches()
        while (len(branches) > 0): # Keep performing until unique rearrangement.
            
            bra  = choice(branches)
            branches.remove(bra)
            enum = topol.iterTypeForBranch(bra,type) # Iterate over each.
        
            for en in enum:
                
                # Get metadata.
                typ = en.getType()
                t   = en.toTree()
                new = t.getNewick()
    
                # See if already been found.
                inlandscape = self.findTreeTopologyByStructure(t.getStructure())
                if (inlandscape != None):
                    # Is already in landscape; has connection to tree?
                    if (inlandscape != i and self.graph.has_node(inlandscape)):
                        if not self.graph.has_edge(inlandscape,i):
                            self.graph.add_edge(inlandscape,i)
                            self.getEdge(inlandscape,i)['weight'] = \
                                self.defaultWeight
                    continue
    
                # See if tree violating existing locks.
                if isvi:
                    en = en.toTopology()
                    if self._isViolating(en): continue
    
                # Score.
                scr = parsimony(new,p)
                t.score  = (None,scr)
                t.origin = typ
                
                # Add to landscape.
                j = self._newNode(t,score=True)
                self.graph.add_edge(i,j)
                self.getEdge(i,j)['weight'] = self.defaultWeight
                return j
        
        # Set explored to True.
        node['explored'] = True 
        
        return None        

    def exploreTree(self,i,type=TYPE_SPR):
        
        ''' Get all neighbors to a tree named i in the landscape using a
        respective rearrangement operator as defined in the rearrangement
        module. Rearrangement type is provided as a rearrangement module type
        definition of form, for example, TYPE_SPR, TYPE_NNI, etc. By default,
        this is TYPE_SPR.
        
        :param i: a tree index
        :param type: the type of rearrangement (e.g., TYPE_SPR, TYPE_NNI)
        :return: a list of neighbors as tree names (usually integers)
        
        '''
        
        # Get parsimony profiles.
        p = self.parsimony_profiles
        
        # Check node.
        node  = self.getNode(i)
        if (node['explored']): return list()
        tre   = node['tree']
        topol = tre.toTopology()
        new   = topol.toNewick()
        isvi  = self._isViolating(topol)
        
        # Check for locks.
        if hasattr(self,'locks'):
            for lock in self.locks:
                hasthis = topol.getBranchFromBipartition(lock)
                if (hasthis): topol.lockBranch(hasthis)
        
        # Perform full-enumeration exploration (1 move).
        neighbors = list()        
        enum = topol.allType(type)
        
        for en in enum:
            
            # Get metadata.
            typ = en.getType()
            t   = en.toTree()
            new = t.getNewick()

            # See if already been found.
            inlandscape = self.findTreeTopologyByStructure(t.getStructure())
            if (inlandscape != None):
                # Is in landscape; has connection to tree?
                if ((inlandscape != i) and self.graph.has_node(inlandscape)): 
                    if not self.graph.has_edge(inlandscape,i):
                        self.graph.add_edge(inlandscape,i)
                        self.getEdge(inlandscape,i)['weight'] = \
                            self.defaultWeight
                continue

            # See if tree violating existing locks.
            if isvi:
                en = en.toTopology()
                if self._isViolating(en): continue

            # Score.
            scr = parsimony(new,p)
            t.score  = (None,scr)
            t.origin = typ
            
            # Add to landscape.
            j = self._newNode(t,score=True)
            neighbors.append(j)
            self.graph.add_edge(i,j)
            self.getEdge(i,j)['weight'] = self.defaultWeight
        
        # Set explored to True.
        node['explored'] = True 
        
        return neighbors

    # Lock Management
        
    def getLocks(self):

        ''' Get all restrictions (locks which cannot be violated on splits).
        
        :return: a list of :class:`.tree.bipartition` objects
        
        '''
        
        return self.locks 
    
    def toggleLock(self,lock):
        
        ''' Add a biparition to the list of locked bipartitions if not 
        present; otherwise, remove it. Return status of lock.
        
        :param lock: a bipartition that cannot be violated
        :type lock: a :class:`.tree.bipartition` object
        :return: a boolean (on or off)
        
        '''
        
        toggle = False
        if lock in self.getLocks():
            self.locks.remove(lock)
        else:
            self.locks.append(lock)
            toggle = True
        return toggle
    
    def lockBranchFoundInTree(self,tr,br):
        
        ''' Given a tree node and a branch object, add a given 
        bipartition to the bipartition lock list. Returns 
        bipartition if locked.
        
        :param tr: a tree name
        :param br: a branch in that tree
        :type br: a :class:`.base.treeBranch` object
        :return: a :class:`.tree.bipartition` object that has been locked or None
        
        '''

        # Check for presence of tree.
        g = self.graph.node
        if (not tr in g): return None
        
        # Acquire metadata for tree.
        tr = self.getTree(tr)
        topl = tr.toTopology()
        
        # Get bipartition for branch.
        bipa  = tree.bipartition(topl,br)
        
        # Toggle the lock.
        self.toggleLock(bipa)
        return bipa
    
    def getBipartitionFoundInTreeByIndex(self,tr,brind,topol=None):
        
        ''' Given a tree node and a branch index, return the associated
        bipartition.
        
        :param tr: a tree name
        :param brind: a branch index in that tree, a la post-order traversal
        :return: a :class:`.tree.bipartition` object
        
        '''
        
        # Check for presence of tree.
        if (not topol) and (not tr in self.graph.node): return None
        
        # Acquire metadata for tree.
        if not topol:
            tr = self.graph.node[tr]['tree']
            topl = tr.toTopology()
        
        # Get bipartition for branch.
        nodes = postOrderTraversal(topl.getRoot())
        if not brind in xrange(len(nodes)): return None
        bnode = nodes[brind]
        bipa  = tree.bipartition(topl,bnode.parent)
        
        # Return the bipartition.
        return bipa
            
    def lockBranchFoundInTreeByIndex(self,tr,brind):
        
        ''' Given a tree node and a branch index, add an associated 
        bipartition to the bipartition lock list. Returns the
        bipartition if locked.
        
        :param tr: a tree name
        :param brind: a branch index
        :return: a :class:`.tree.bipartition` object if it has been locked
        
        '''

        # Check for presence of tree.
        if (not tr in self.graph.node): return None
        
        # Acquire metadata for tree.
        tree = self.graph.node[tr]['tree']
        topl = tree.toTopology()
        
        # Get bipartition for branch.
        bipa = self.getBipartitionFoundInTreeByIndex(tr,brind,topl)
        
        # Toggle the lock.
        self.toggleLock(bipa)
        return bipa

    def _isViolating(self,topo):
        
        ''' PRIVATE: Given a topology, determine if it is
        violating existing locks. '''
        
        if len(self.locks) > 0:
            for lock in self.locks:
                hasit = topo.getBranchFromBipartition(lock)
                if (not hasit): return True
        return False
    
    def isViolating(self,i):
        
        ''' Determine if a tree is violating any locks intrinsic 
        to the landscape. Will also return False if the tree (name) is not
        present in the landscape.
        
        :param i: a tree name (usually an integer)
        
        '''
    
        if i in self.graph.node.keys():
            node = self.getNode(i)
            topo = node['tree'].toTopology()
            return self._isViolating(topo)
        return False

    # Search Methods

    def __getitem__(self, i):
        
        return self.graph.node[i]

    def indexOf(self, tr):
        
        ''' Acquire the index/name in this landscape of a tree object. Returns
        -1 if not found. Warning: naively performs a sequential search.
        
        :param tr: a tree
        :type tr: a :class:`.tree.tree` object
        :return: a tree name (usually an integer index) or -1 if not found
        
        '''
        
        for t in self.graph.node:
            if self.getTree(t) == tr: return t
        return -1

    def findTree(self,newick):
        
        ''' Find a tree by Newick string, taking into account branch lengths.
        Returns the index of this tree in the landscape. Warning: naively
        performs a sequential search.
        
        :param newick: a Newick string
        :type newick: a string
        :return: a tree name (usually an integer index) or None if not found
        
        '''
        
        for t in self.graph.node:
            if self.getTree(t).newick == newick: return t
        return None
    
    def findTreeTopology(self,newick):
        
        ''' Find a tree by topology, not taking into account branch lengths.
        
        :param newick: a Newick string
        :type newick: a string
        :return: a tree name (usually an integer index) or None if not found
        
        '''
        
        s = newickParser(newick).parse()
        removeBranchLengths(s)
        s = str(s) + ';'
        return self.findTreeTopologyByStructure(s)
    
    def findTreeTopologyByStructure(self,struct):
        
        ''' Find a tree by topology, not taking into account branch lengths,
        given the topology.
        
        :param struct: a Newick string without branch lengths (a "structure")
        :type struct: a string
        :return: a tree name (usually an integer index) or None if not found
        
        '''
        
        query = (struct in self.newickSearchDict)
        if query is True:
            index = self.newickSearchDict[struct]
            return index
        return None
        
    def getBestImprovement(self,i):
        
        ''' For a tree in the landscape, investigate neighbors to find 
        a tree that leads to the best improvement of fitness function score
        on the basis of likelihood.
        
        :param i: a tree name (usually an integer)
        :return: a tree name (usually an integer) or None if no better tree
        
        '''
        
        nodes = self.getNode
        tree  = self.getTree
        ml    = lambda d: tree(d).getScore()[0]
        if (not i in self.getNodeNames()):
            raise LookupError('No tree by that name in landscape.' % (i))
        node  = nodes(i)
        near  = self.graph.neighbors(i)
        if (len(near) == 0): return None
        best  = max(near,key=ml)
        if (best != None and ml(best) > ml(i)):
            return best
        else: return None
        
    def getPathOfBestImprovement(self,i):
        
        ''' For a tree in the landscape, investigate neighbors iteratively until
        a best path of score improvement is found on basis of likelihood.
        
        :param i: a tree name (usually an integer)
        :return: a list of tree names (usually integers)
        
        '''
        
        path  = list()
        node = self.getNode(i)
        impr = self.getBestImprovement(i)
        curs = node
        if impr == None: return path
        
        while (impr != None):
            path.append(impr)
            impr = self.getBestImprovement(impr)
        
        return path
    
    def getAllPathsOfBestImprovement(self):
        
        ''' Return all paths of best improvement as a dictionary.
        
        :return: a dictionary of tree name to paths (lists of tree names)
        
        '''
        
        paths = dict()
        for node in self.graph.node:
            paths[node] = self.getPathOfBestImprovement(node)
        return paths    
    
    def iterAllPathsOfBestImprovement(self):
        
        ''' Return an iterator for all paths of best improvement.
        
        :return: a generator of paths (lists of tree names)
        
        '''
        
        for node in self.graph.node:
            yield self.getPathOfBestImprovement(node)
    
    def isLocalOptimum(self,i):
        
        ''' Determine if a tree is a local optimum. This means it has the
        following properties: 
        
          (1) Possesses a likelihood score.
          (2) Local neighborhood completely enumerated (and scored).
          (3) None of its neighbors is a better improvement.
        
        :return: a boolean
        
        '''
        
        nodes = self.graph.node
        if (self.getTree(i).score[0] == None):
            return False
        elif (not self.getNode(i)['explored']):
            return False
        p = self.getPathOfBestImprovement(i)
        if (len(p) == 0):
            near = self.graph.neighbors(i)
            if (len(near) == 0): return False
            for node in [nodes[x] for x in near]:
                if node['tree'].score[0] == None and \
                   'failed' in node and not node['failed']:
                    return False
            return True
        return False
    
    def getLocalOptima(self):
        
        ''' Get all trees in the landscape that can be labelled as a local
        optimum.
        
        :return: a list of tree names (usually integers)
        
        '''
        
        opt = []
        for node in self.graph.node:
            if (self.isLocalOptimum(node)): opt.append(node)
        return opt
    
    def getGlobalOptimum(self):
        
        ''' Get the global optimum of the current space.
        
        :return: a tree name (usually an integer)
        
        '''
        
        # TODO: store global optimum as trees are added?
        optima = self.getLocalOptima()
        return max(optima,key=lambda d: self.getTree(d).score[0])

    # Output Methods
    
    def __str__(self): return self._dump()
    
    def _dump(self):
        
        trees = []
        for t in self.getNodeNames(): trees.append(self.getVertex(t))
        if type(self.alignment) == alignment.phylipFriendlyAlignment:
            return '\n'.join([t.getProperNewick() for t in trees])
        else: return '\n'.join([t.getNewick() for t in trees])

    def _toTreeSet(self, func):
        
        treeset = treeSet()
        for t in self.getNodeNames():
            treeset.addTree(func(t))
        return treeset

    def toProperNewickTreeSet(self):
        
        ''' Convert this landscape into an unorganized set of trees where taxa
        names are transformed to their original form (i.e. not transformed to a
        state friendly for the Phylip format).
        
        :return: a :class:`.tree.treeSet` object
        
        '''
        
        return self._toTreeSet(
            lambda t: tree(self.getVertex(t).getProperNewick()))
    
    def toTreeSet(self):

        ''' Convert this landscape into an unorganized set of trees.
        
        :return: a :class:`.tree.treeSet` object
        
        '''
        
        return self._toTreeSet(lambda t: self.getTree(t))

# Comprising Vertices of Landscapes


class vertex(object):
    
    ''' Encapsulate a single vertex in the landscape and add convenient
    functionality to alias parent landscape functions. '''
    
    def __init__(self, obj, ls):
        
        ''' Initialize this vertex. '''
        
        self.id  = obj['index']
        self.obj = obj
        self.ls = ls
        self.ali = ls.alignment
    
    def getIndex(self):
        
        ''' Get the index of this tree in the space.
        
        :return: a tree name (usually an integer)
        
        '''
        
        return self.id
    
    def getDict(self):
        
        ''' Get the dictionary object (key-value pairs) associated with this
        tree as it is in the NetworkX graph. 
        
        :return: a dictionary
        
        '''
        
        return self.obj
    
    def getObject(self):
        
        ''' Get the dictionary object (key-value pairs) associated with this
        tree as it is in the NetworkX graph. 

        :return: a dictionary

        '''
        
        return self.getDict()
    
    def getTree(self):
        
        ''' Get the tree object associated with this tree. 
        
        :return: a :class:`.tree.tree` object
        
        '''
        
        return self.obj['tree']
    
    def getNewick(self):
        
        ''' Get the Newick string of this tree. 
        
        :return: a string
        
        '''
        
        return self.getTree().getNewick()
    
    def getScore(self):
        
        ''' Get (any) score(s) associated with this tree. 
        
        :return: a tuple of floating point values (scores)
        
        '''
        
        return self.getTree().score
    
    def getOrigin(self):
        
        ''' Get the origin of this tree (how it was acquired).
        
        :return: a string
        
        '''
        
        return self.getTree().origin
    
    def getNeighbors(self):
        
        ''' Get any neighbors to this tree in the landscape. 
        
        :return: a list of tree names (usually integers)
        
        '''
        
        return self.ls.graph.neighbors(self.id)
    
    def getDegree(self):
        
        ''' Get the degree of this tree in the graph. 
        
        :return: an integer
        
        '''
        
        return len(self.getNeighbors())
    
    def isLocalOptimum(self):
        
        ''' Determine if this tree is an optimum. 
        
        :return: a boolean
        
        '''
        
        return self.ls.isLocalOptimum(self.id)
    
    def isExplored(self):
        
        ''' See if this tree has had all possible rearrangements performed. 
        
        :return: a boolean
        
        '''
        
        return self.obj['explored']
    
    def isFailed(self):
        
        ''' Determine if any errors are associated with this node.
        
        :return: a boolean
        
        '''
        
        return ('failed' in self.obj and self.obj['failed'])
    
    def setExplored(self, exp):
        
        ''' Override the "explored" flag of this node in the landscape.
        
        :param exp: a boolean
        
        '''
        
        self.obj['explored'] = exp
    
    def approximatePossibleNumNeighbors(self):
        
        ''' Approximate the possible number of neighbors to this vertex
        by considering the type of tree rearrangement operator. Returns
        LS_NOT_DEFINED if the operator is not known yet.
        
        :return: an integer
        
        '''

        n = self.ls.getNumberTaxa()
        if self.ls.operator == 'SPR': return 4*(n-3)*(n-2)
        return LS_NOT_DEFINED
    
    def scoreLikelihood(self):
    
        ''' Acquire the log-likelihood for this vertex.
        
        :return: the log-likelihood score
        
        '''
        
        if self.getScore()[0] == None:
            l = ll(self.getTree(),self.ali)
            self.obj['tree'].score = (l,self.getScore()[1])
            return l
        else: return self.getScore()[0]
    
    def getBestImprovement(self):
        
        ''' Alias function for function of same name in parent
        landscape. '''
        
        return self.ls.getBestImprovement(self.id)
    
    def getPathOfBestImprovement(self):
    
        ''' Alias function for function of same name in parent
        landscape. '''
        
        return self.ls.getPathOfBestImprovement(self.id)
    
    def isBestImprovement(self):
        
        ''' Check to see if this vertex is a best move for another node.
        
        :return: a boolean
        
        '''
        
        for neighbor in self.getNeighbors():
            if self.ls.getTree(neighbor).score[0] < self.getTree().score[0]:
                if self.ls.getBestImprovement(neighbor) == self.id:
                    return True
        return False
    
    def isViolating(self):
        
        ''' Alias function for function of same name in parent
        landscape. '''
        
        return self.ls.isViolating(self.id)

    def getProperNewick(self):
        
        ''' Get the proper Newick string for a tree. 
        :returns: A string.'''
        
        if self.ali == None: return self.getNewick()
        return self.ali.reinterpretNewick(self.getNewick())
    
    def iterBipartitions(self):
        
        ''' Return a generator to iterate over all bipartitions for this vertex.
        '''
    
        # Get tree object.
        tr = self.getTree()
        
        # Get corresponding topology.
        topo  = tr.toTopology()
        nodes = postOrderTraversal(topo.getRoot())
        
        # Get branches + bipartitions of the topology.
        bp = []
        br = [x.parent for x in nodes]
        for b in br:
            if b:
                bi = tree.bipartition(topo,b)
                yield bi
            else: yield None
    
    def getBipartitions(self):
    
        ''' Get all bipartitions for this vertex.
        
        :return: a list of :class:`.tree.bipartition` objects
        
        '''
        
        # Return the bipartitions.
        return [bipart for bipart in self.iterBipartitions()]
    
    def getBipartitionScores(self):
        
        ''' Get all corresponding bipartition vectors of SPR scores. '''
        
        return [x.getSPRScores() for x in self.getBipartitions()]
    
    def getNeighborsOfBipartition(self, bi):
        
        ''' Get corresponding neighbors of a bipartition in this vertex's tree.
        '''
        
        l     = self.ls
        re    = bi.getSPRRearrangements()
        lsids = [l.findTreeTopology(x.toTree().getStructure()) for x in re]
        return lsids
    
    def getNeighborsOfBranch(self, br):
        
        ''' Get corresponding neighbors of a branch in this vertex's tree. '''
        
        tr   = self.getTree()
        topl = tr.toTopology()
        bipa = tree.bipartition(topl,br)
        return self.getNeighborsOfBipartition(bipa)

    
