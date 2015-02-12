''' Phylogenetic tree structure encapsulation; allow rearrangement of said structure. Tree rearrangements inducing other topologies include Nearest Neighbor Interchange (NNI), Subtree Pruning and Regrafting (SPR), and Tree Bisection and Reconstruction (TBR). Each of these describe a transfer of one node in phylogenetic trees from one parent of a tree to a new parent. Respectively, these operators describe transformations that are subsets of those possible by the successive operator. For example, an NNI operator can perform transformations that are a subset of the transformations possible by the SPR operator. '''

# Date:   Feb 10 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import newick, tree, base

# Exception Handling

class RearrangementError(Exception):
    def __init__(self,val): self.value = val
    def __str__(self): return repr(self.value)

# Movement Type Definitions

TYPE_SPR, TYPE_NNI, TYPE_TBR = 1, 2, 3

# Simulate deep copying.

def dup(topo,where=None):
    if where: new = topology(toLeaf=where)
    else: new = topology()
    new.fromNewick(topo.toNewick())
    if where:
        fake = tree.bipartition(topo,new.fakebranch)
        # Retain locks.
        for lock in new.locked:
            if lock != fake: new.locked.append(lock)      
    return new

# Rearrangement structure.

class rearrangement:
    
    ''' Encapsulates a single rearrangement move
    of type SPR, NNI, ... '''
    
    def __init__(self,struct,type,targ,dest):
        
        ''' Initialize by providing a pointer to
        a base topology, a target branch 
        to be moved, and its destination. '''
        
        self.topol       = struct
        self.type        = type
        self.target      = targ
        self.destination = dest    
    
    def getType(self):
        
        ''' Get the type of movement. '''
        
        if (self.isSPR()):   type = 'SPR'
        elif (self.isNNI()): type = 'NNI'
        elif (self.isTBR()): type = 'TBR'
        else:                type = 'UNKNOWN'
        return type
    
    def isNNI(self): return (self.type == TYPE_NNI)
    def isSPR(self): return (self.type == TYPE_SPR)
    def isTBR(self): return (self.type == TYPE_TBR)
    
    def toTopology(self):
        
        ''' Commit the actual move and return the
        topology. '''
        
        topo = self.topol.move(self.target,self.destination)
        topo.rerootToLeaf()
        return topo
    
    def toNewick(self):
        
        ''' Commit the move but do not create a 
        new structure. Only retrieve resultant
        Newick string; will be more efficient.'''
        
        return self.topol.move(self.target,self.destination,
                               returnStruct=False)
    
    def toTree(self):
        
        ''' Commit the move and transform to tree object. '''
        
        out = self.toTopology().toTree()
        out.origin = self.getType()
        return out
    
    def doMove(self): return self.toTopology()
    
    def __str__(self):
        
        t = self.getType()
        return '<%s> move [%s] to [%s]' % (t,self.target,self.destination)

# Tree structure.

class topology(base.treeStructure):
    
    ''' Encapsulate a tree topology, wrapping the newick 
    tree structure. Is immutable. '''
    
    def __init__(self,t=None,rerootToLeaf=True,toLeaf=None):

        ''' Initialize structure with a top-level internal 
        node OR nothing. '''
        
        self.root        = t
        self.forbidden  = {}
        self.branches   = []
        self.locked     = []
        self.fakebranch = None
        self.rerootFlag = rerootToLeaf
        self.rerootLoc  = toLeaf
        
        if t != None: 
            if rerootToLeaf: self.rerootToLeaf(toLeaf)
            self._getAllBranches()
            self._getForbiddenStates()
            self._clearInteriorNodeNames()
            if rerootToLeaf: self._lockLeafBranch()       
       
    def _getAllBranches(self):
        
        ''' PRIVATE: Populate branch list. '''
        
        for br in self.root.children:
            self.branches.append(br)
            self.branches.extend(newick.getAllBranches(br))
        self.partitions = [[x for x in self.branches]]  
       
    def _getForbiddenStates(self):
        
        ''' PRIVATE: Construct forbidden state lists for 
        each branch. '''
        
        # Set up list for each branch.
        for br in self.branches:
            # Anything inside the subtree for that branch
            # is forbidden to be moved to by this branch.
            self.forbidden[br] = [_ for _ in newick.getAllBranches(br)]
            # Any siblings are also forbidden 
            # (give the same tree).
            siblings = br.parent.children
            siblings = [_ for _ in siblings if _!=br]
            self.forbidden[br].extend(siblings)
            # And the parent is forbidden.
            self.forbidden[br].append(br.parent.parent)
        
    def _clearInteriorNodeNames(self):
        
        ''' PRIVATE: Clear all interior node names 
        but store them in case are needed. '''
        
        # Make sure interior nodes are unnamed.
        nodes = self.getAllNodes()
        for n in nodes:
            if n != self.root and len(n.children) != 0 and n.label != '':
                n._label = n.label
                n.label = '' 

    def _lockLeafBranch(self):
        
        ''' PRIVATE: Lock the branch hanging 
        off the root with a leaf. '''
        
        # Get this branch.
        lbranch = None
        for f in self.root.children:
            if (len(f.child.children) == 0):
                lbranch = f
        if (lbranch == None):
            raise RearrangementError(
                'Could not locate leaf branch off root to lock.')
        self.lockBranch(lbranch)

    def _getPartition(self,b):
        
        ''' PRIVATE: Acquire the list of branches: those this particular
        branch could move to without violation of any locks 
        and itself. '''
        
        possible = [x for x in self.branches]
        for lo in self.locked:
            l,r = lo.getBranchListRepresentation()
            if b in l: o = r
            else:      o = l
            for i in o:
                if i in possible: possible.remove(i)
        return possible

    def _flipOperation(self,opname,br):
        
        ''' PRIVATE: Try a particular rearrangement operation on a flipped
        version of the tree. '''

        # Check what bipartition is represented.
        bipart = tree.bipartition(self,br)
        # Get a leaf that can be rerooted to.
        newleaf = None
        leaves  = base.treeStructure.leaves(br.child)
        curleaf = sorted(
            self.getAllLeaves(),key=lambda d: d.label)[0]
        for leaf in leaves:
            if leaf != curleaf:
                newleaf = leaf
                break
        if (newleaf == None): return []
        # Make flipped tree structure.
        n = dup(self,newleaf)
        # Get corresponding branch.
        b = n.getBranchFromBipartition(bipart)
        # Do possible moves from the other way.
        op = getattr(n,opname)
        return op(b,flip=False)

    def rerootToLeaf(self,toleaf=None):
        
        ''' PRIVATE: Reroots the given tree structure such
        that it is rooted nearest the lowest-order leaf. '''
        
        # Determine lowest-order leaf.
        if not toleaf:
            toleaf = sorted(self.getAllLeaves(),
                            key=lambda d: d.label)[0]
        else:
            # Find it in current topology.
            found = False
            for leaf in self.getAllLeaves():
                if leaf.label == toleaf.label:
                    toleaf = leaf
                    found = True
                    break
            if (not found): raise RearrangementError(
                'While rerooting, could not find leaf to reroot to.')

        # Acquire branch connecting to leaf + adjacent node.
        branch  = toleaf.parent
        closest = branch.parent
        
        # Error checking: either of these two do not exist?
        if not (branch):  raise RearrangementError(
            'While rerooting, leaf has no connecting branch.')
        if not (closest): raise RearrangementError(
            'While rerooting, connecting branch of leaf has no parent.')

        # Flip all of the directionality to this node.
        if closest != self.root: 
            t = newick.invertAlongPathToNode(closest,self.root)
            if not t: raise RearrangementError(
                'While rerooting, directions could not be reversed.')
        
        # Strip the lowest node off of the tree.
        closest.children.remove(branch)
        
        # Create a fake branch.
        fakebr = newick.branch(closest,0.0)
        closest.parent = fakebr
        
        # Create a new node to act as the new root.
        root = newick.node(children=[branch,fakebr])
        branch.parent = root
        fakebr.parent = root
        
        # Reassign top-level node.
        self.root = root       
        
        # Remove unary internal nodes.
        newick.removeUnaryInternalNodes(self.root)
        
        # Ensure fake branch is recognized.
        for f in root.children:
            if (f.branch_length == 0 and len(
                f.child.children) > 0):
                fakebr = f
        self.fakebranch = fakebr      
    
    def getBranches(self): return self.branches
    def getLeaves(self): return self.getAllLeaves()
    
    def getBipartitions(self):
        
        ''' Get all bipartitions. '''
        
        bili = []
        br = self.branches
        for b in br:
            bi = tree.bipartition(self,b)
            if not bi in bili: bili.append(bi)
        return bili
    
    def getStrBipartitionFromBranch(self,br):
        
        ''' Given a branch, return corresponding 
        bipartition. '''
        
        right  = base.treeStructure.leaves(br.child)
        others = self.getAllLeaves()
        r      = [x.label for x in right]
        l      = [x.label for x in others 
                  if not x.label in r]
        return l,r
    
    def getBranchFromStrBipartition(self,bip):
        
        ''' Given a bipartition of taxa, return a 
        branch that creates that partition of tree
        taxa. '''
        
        l,r = sorted(bip[0]),sorted(bip[1])
        
        for br in self.branches:
            right  = base.treeStructure.leaves(br.child)
            names  = sorted([x.label for x in right])
            if (names == l) or (names == r):
                return br
        return None
    
    def getBranchFromBipartition(self,bip):
        
        ''' Given a bipartition object, return
        a branch that creates that partition of
        taxa. '''
        
        return self.getBranchFromStrBipartition(
            bip.getStringRepresentation())

    def lockBranch(self,branch):
        
        ''' Given a branch, lock it such that no 
        transitions can ever occur across it. '''
        
        # Transform to bipartition object.
        bipart = tree.bipartition(self,branch)
        
        # Has already been locked?
        if bipart in self.locked: return True
            
        # Is it even in this topology?
        isintopol = (branch in self.branches)
        if (not isintopol): return False
        self.locked.append(bipart)        
        return True

    def move(self,branch,destination,returnStruct=True):
        
        ''' Move a branch and attach to a destination 
        branch. Return new structure, or return merely 
        the resultant Newick string. '''
        
        # Cannot move to these.
        forbidden = self.forbidden[branch]
        if destination in forbidden:
            raise RearrangementError(
                'Cannot move into subtree or to sibling.')
        
        # Check partitions.
        partition = self._getPartition(branch)
        if destination not in partition:
            raise RearrangementError(
                'Cannot move outside of locked partition.')
        
        # Get parents.
        s_parent = branch.parent       # Source parent. 
        t_parent = destination.parent  # Target parent.
        
        # If parents were found, do transform.
        if s_parent != None and t_parent != None:
            
            # Remove the branch.
            s_parent.children.remove(branch)
            
            # Create new node.
            node = newick.node('',[branch])
            
            # Break destination branch in half, attach.
            half  = destination.branch_length/2.0
            outer = newick.branch(destination.child,half,node)
            inner = newick.branch(node,half,t_parent)
            t_parent.children.remove(destination)
            t_parent.children.append(inner)
            node.children.append(outer)
            destination.child.parent = outer
            node.parent = inner
            
            # Check degree of source parent; combine edges if necessary.
            if len(s_parent.children) == 1:
                a           = s_parent.children[0]
                b           = s_parent.parent
                end         = a.child
                start       = b.parent
                comb        = a.branch_length + b.branch_length
                fakebr      = newick.branch(end,comb,start)
                end.parent  = fakebr
                start.children.remove(b)
                start.children.append(fakebr)
            
            if returnStruct:
                # Immutable so recreate new structure.
                result = dup(self)
            else: result = self.toNewick()
            
            # Undo changes.
            s_parent.children.append(branch)
            if len(s_parent.children) == 2:
                end.parent  = a
                start.children.append(b)
                start.children.remove(fakebr)
            t_parent.children.append(destination)
            t_parent.children.remove(inner)
            destination.child.parent = destination

            return result
        
        else: return None
        
    def SPR(self,branch,destination):
        
        ''' Perform an SPR move of a branch to a destination 
        branch, creating a new node there. Returns a 
        rearrangement structure (not the actual new structure) 
        that can then be polled for the actual move; this 
        is in order to save memory. '''
        
        return rearrangement(self,TYPE_SPR,branch,destination)
    
    def NNI(self,branch,destination):
        
        ''' Perform an NNI move of a branch to a destination, 
        only if that destination branch is a parent's parent 
        or a parent's sibling. Returns a rearrangement structure
        (not the actual new structure) that can then be polled
        for the actual move; this is in order to save memory. '''
        
        if newick.isSibling(branch.parent,destination) or \
           (branch.parent and destination == 
            branch.parent.parent):
            return rearrangement(self,TYPE_NNI,branch,destination)
        else: return None
    
    def iterSPRForBranch(self,br,flip=True):    
        
        ''' Consider all valid SPR moves for a given branch 
        in the topology and yield all possible rearrangements
        as a generator. '''    
    
        # Go through possible moves.
        partition = self._getPartition(br)
        if (partition == None):
            raise RearrangementError(
                'Branch not in any partition of topology.')
        forbidden = self.forbidden[br]
        possible  = [x for x in partition if not
                     x in forbidden and not
                     x == br]
        
        # Pass rearrangement structure as yielded object.
        for dest in possible:
            move = self.SPR(br,dest)
            if (move): yield move

        # Flip tree around and try other way.
        if flip:
            ite = self._flipOperation('iterSPRForBranch',br) 
            for it in ite: yield it
    
    def allSPRForBranch(self,br,flip=True):
        
        ''' Consider all valid SPR moves for a given branch 
        in the topology and return all possible rearrangements. '''
        
        return [x for x in self.iterSPRForBranch(br,flip)]
    
    def allSPR(self):
        
        ''' Consider all valid SPR moves for a given topology
        and return all possible rearrangements. '''

        # Output list of structures.
        li = []        
        
        for branch in self.branches:
            li.extend(self.allSPRForBranch(branch))
        
        # Return the list.
        return li

    def iterNNIForBranch(self,br,flip=True):

        ''' Consider all valid NNI moves for a given branch
        in the topology and and yield all possible rearrangements
        as a generator. '''
        
        # Go through possible moves.
        partition = self._getPartition(br)
        if (partition == None):
            raise RearrangementError(
                'Branch not in any partition of topology.')
        possible  = [br.parent.parent] + br.parent.children
        forbidden = self.forbidden[br]
        possible  = [x for x in partition if x in possible and 
                     not x in forbidden and not x == None]
        
        # Pass rearrangement structure as yielded object.
        for dest in possible:
            move = self.NNI(br,dest)
            if (move): yield move

        # Flip tree around and try other way.
        if flip:
            ite = self._flipOperation('iterNNIForBranch',br) 
            for it in ite: yield it        

    def allNNIForBranch(self,br,flip=True):
        
        ''' Consider all valid NNI moves for a given branch 
        in the topology and return all possible rearrangements. '''
        
        return [x for x in self.iterNNIForBranch(br,flip)]

    def allNNI(self):
        
        ''' Consider all valid NNI moves for a given topology
        and return all possible rearrangements. '''

        # Output list of structures.
        li,nni = [],self.allNNIForBranch        
        for branch in self.branches: li.extend(nni(branch))
        return li
    
    def allType(self,type=TYPE_SPR):
    
        ''' Consider all valid moves of a given rearrangement
        operator for a given topology. Uses a given rearrangement 
        operator type defined in this module. For example, calling this 
        function by providing TYPE_NNI as the type will iterate over all 
        NNI operations. By default, the type is TYPE_SPR. '''
        
        if (type == TYPE_SPR):   return self.allSPR()
        elif (type == TYPE_NNI): return self.allNNI()
        else: raise RearrangementError('No rearrangement type of that form is defined.')
    
    def iterTypeForBranch(self,br,type=TYPE_SPR,flip=True):
        
        ''' Iterate over all possible rearrangements for a
        branch using a given rearrangement operator type defined
        in this module. For example, calling this function by providing
        TYPE_NNI as the type will iterate over all NNI operations. By
        default, the type is TYPE_SPR. '''
        
        if (type == TYPE_SPR):   return self.iterSPRForBranch(br,flip)
        elif (type == TYPE_NNI): return self.iterNNIForBranch(br,flip)
        else: raise RearrangementError('No rearrangement type of that form is defined.')
    
    def fromNewick(self,newickstr):
        
        ''' Alias for parse(). '''
        
        return self.parse(newickstr)
    
    def parse(self,newickstr):
        
        ''' Parse a newick string and assign the tree to this
        object. Cannot already be initialized with a tree. '''
    
        if self.root != None:
            raise RearrangementError(
                'Structure already initialized.')
        p = newick.newickParser(newickstr)
        self.root  = p.parse()
        self.orig = newickstr
        if self.rerootFlag: self.rerootToLeaf(self.rerootLoc)
        self._getAllBranches()
        self._getForbiddenStates()
        self._clearInteriorNodeNames()
        if self.rerootFlag: self._lockLeafBranch()
        
    def toNewick(self):
        
        ''' Return the newick string of the tree. '''
        
        return self.__str__()
    
    def toUnrootedNewick(self):
        
        ''' Return the newick string of the tree as an
        unrooted topology with a multifurcating top-level node. '''
        
        r,a,b = self.root,self.root.children[0],self.root.children[1]
        r.children.remove(a)
        b.child.children.append(a)
        a.parent = b.child
        newic = str(b.child) + ';'
        r.children.append(a)
        b.child.children.remove(a)
        a.parent = self.root
        return newic
        
    def toTree(self):
        
        ''' Return the tree object for this topology. '''
        
        return tree.tree(self.toNewick())
        
    def toUnrootedTree(self):
        
        ''' Return the tree object of the unrooted 
        version of this topology. '''
        
        return tree.tree(self.toUnrootedNewick())
        
    def __str__(self):
        
        ''' Return the newick string of the tree. '''
        
        return str(self.root) + ';'
