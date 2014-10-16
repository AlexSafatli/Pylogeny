''' Newick string parsing and object interaction. A Newick string can represent a phylogenetic tree. '''

# Date:   Feb 16 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import rearrangement
from random import shuffle
from math import factorial as fact

# Exception Handling

class ParsingError(Exception):
    def __init__(self,val): self.value = val
    def __str__(self): return repr(self.value)

# Possible Number of Trees

numberRootedTrees   = lambda t: numberUnrootedTrees(t+1)
numberUnrootedTrees = lambda t: (fact(2*(t-1)-3))/((2**(t-3))*fact(t-3))

# Class Definitions

class tree(object):
    
    ''' Defines a single phylogenetic tree by newick string;
    can possess other metadata. '''
    
    def __init__(self,newi='',check=False):
        
        self.name   = ''
        self.score  = None
        self.origin = None
        self.newick = newi
        if (check): self._checkNewick(newi)        
        else: self._setNewick(newi)

    # Internals
        
    def __eq__(self,o): return (o.struct == self.struct)
    def __ne__(self,o): return not (self.__eq__(o))
    def __str__(self): return self.newick
    
    def _checkNewick(self,newi):
        
        ''' PRIVATE: Run the Newick string through a 
        parse pass and reroot to lowest-order leaf in
        order to ensure a consistent Newick string. '''
        
        newi        = self._cleanNewick(newi)
        prsd        = self._getParsed(newi)
        topo        = rearrangement.topology(prsd)
        self.newick = topo.toNewick()
        self.struct = self._getStructure(prsd)
    
    def _getParsed(self,newi):
        
        ''' PRIVATE: Acquires a topology structure for
        input Newick string. '''
        
        return parser(newi).parse()
    
    def _getStructure(self,prsd=None):
        
        ''' PRIVATE: Acquires a newick string without any
        defined branch lengths. '''
        
        if prsd: p = prsd
        else: p = self._getParsed(self.newick)    
        removeBranchLengths(p)
        return str(p) + ';'
        
    def _cleanNewick(self,newi):
        
        ''' PRIVATE: Cleans input newick string such that 
        a semicolon is added if not present, etc. '''
        
        return newi.strip('\n').strip(';') + ';'   
    
    def _setNewick(self,n):
        
        ''' PRIVATE: Set Newick string with little
        intervention. '''
        
        self.newick = n
        self.struct = self._getStructure()    
    
    # Mutators
    
    def setName(self,n):    self.name = n
    def setOrigin(self,o):  self.origin = o    
    def setScore(self,s):   self.score = s
    
    # Accessors
    
    def getName(self):      return self.name
    def getScore(self):     return self.score
    def getOrigin(self):    return self.origin
    def getNewick(self):    return self.__str__()
    def getStructure(self): return self.struct  
    
    def getSimpleNewick(self):    

        ''' Return a Newick string with all Taxa name replaced with
        successive integers. '''
    
        t = self.toTopology()
        n = getAllLeaves(t.getRoot())
        for _ in xrange(len(n)): n[_].label = str(_+1)
        return str(t.getRoot()) + ';'
    
    def toTopology(self):
        
        ''' Return a topology instance for this tree. '''
    
        t = rearrangement.topology()
        t.fromNewick(self.newick)
        return t

class node(object):
    
    ''' Newick node. '''
    
    def __init__(self,lbl='',strees=None,parent=None):
        self.label    = lbl
        self.parent   = parent
        if strees:
            self.children = strees
        else: self.children = list()
    def __str__(self):
        sT, sl = self.children, self.label
        if len(self.children) > 0:
            # Perform sorting to ensure consistent naming scheme.
            l = sorted([st for st in sT],key=lambda d: d.child.label)
            return '(%s)%s' % (','.join(map(str,l)),sl)
        else: return self.label

class branch(object):
    
    ''' Newick branch. '''
    
    def __init__(self,chi,l,parent=None,s=None):
        self.parent         = parent
        self.child          = chi
        self.branch_length  = l
        self.branch_support = s
    def __str__(self):
        st, sl = self.child, self.branch_length
        if sl > 0: return '%s:%s' % (str(st),str(sl))
        else:      return '%s'    % (str(st))
        
# Traversal Functions

def assignParents(top):
    
    ''' Should be a one-time use function. Goes through
    and assigns parents to the parsed newick tree structure
    nodes and branches to allow for up-traversal. '''
    
    # Is a node?
    if type(top) != node: return False
    
    # Traverse.
    for br in top.children:
        br.parent       = top
        br.child.parent = br
        assignParents(br.child)
    return True

def removeBranchLengths(top):
    
    ''' Goes through and removes any stored branch lengths. '''

    # Is a node?
    if type(top) != node: return False
    
    # Traverse.
    for br in top.children:
        br.branch_length = 0.0
        removeBranchLengths(br.child)
    return True

def removeUnaryInternalNodes(top):
    
    ''' Goes through and ensures any degree-2 internal 
    nodes are smoothed into a single degree-3 internal 
    node. '''
    
    next_item = None
    if (len(top.children) == 1):
        pa             = top.parent
        ch             = top.children[0]
        if (pa and ch):
            st, end = pa.parent, ch.child
            st.children.remove(pa)
            br = branch(end,pa.branch_length+ch.branch_length,st)
            st.children.append(br)
            end.parent = br
            next_item  = end
    
    if not next_item:
        for b in top.children:
            removeUnaryInternalNodes(b.child)
    else: removeUnaryInternalNodes(next_item)
    
def invertAlongPathToNode(target,top):
    
    ''' DANGEROUS: Reverses all directionality to a given
    node from a top-level node. Intended as a low-level 
    function for rerooting a tree. '''
        
    def invertDirections(n,newparent=None):
        b = n.parent
        if b:
            b.child  = b.parent
            b.parent = n
            if len(n.children) != 0:
                n.children.append(b)
        n.parent = newparent        
        
    # Are nodes?
    if type(target) != node: return False
    elif type(top)  != node: return False
    
    if target == top:
        # End condition.
        invertDirections(target)
        return True
    
    for b in top.children:
        endpoint = b.child
        if invertAlongPathToNode(target,endpoint):
            invertDirections(top,b)
            top.children.remove(b)
            return True
    return False

def isLeaf(n):
    
    ''' Given a node, see if a leaf. '''
    
    return (len(n.children) == 0)

def isInternalNode(n):
    
    ''' Given a node, see if is an internal node. '''
    
    return (len(n.children) > 1)

def shuffleLeaves(top):
    
    ''' DANGEROUS: Given a top-level node, shuffle all 
    leaves in this tree. '''

    leaves = getAllLeaves(top)
    names  = [x.label for x in leaves]
    shuffle(names)
    for l in xrange(len(leaves)):
        leaves[l].label = names[l]

def getAllLeaves(top):
    
    ''' Given a top-level node, find all leaves. '''

    # Is a node? Is a leaf?
    if type(top) != node: return list()
    if isLeaf(top):       return [top]
    
    # Get list.
    li = list()
    
    # Traverse.
    for br in top.children: li.extend(getAllLeaves(br.child))
    return li

def getAllInternalNodes(top):
    
    ''' Given a top-level node, find all internal nodes. '''
    
    # Is a node? Leaf?
    if type(top) != node:   return list()
    if isLeaf(top):         return list()
    
    # Get list.
    li = list()
    
    # Traverse.
    for br in top.children: li.extend(getAllInternalNodes(br.child))
    return li

def getAllNodes(top):
    
    ''' Given a node, traverse all nodes and return 
    as a list in pre-order. '''

    # Is a node?
    if type(top) != node: return list()

    # Make list.
    li = list()
    
    # Traverse.
    li.append(top)
    for br in top.children: li.extend(getAllNodes(br.child))
    return li

def postOrderTraversal(top):
    
    ''' Given a node, traverse all nodes and return as 
    a list in post-order. '''

    # Is a node?
    if (top == None or type(top) != node): return list()
    
    # Make list.
    li = list()
    
    # Traverse.
    for br in top.children: li.extend(postOrderTraversal(br.child))
    li.append(top)
    return li

def getAllBranches(br):
    
    ''' Given a branch, traverse subtree and return 
    comprising branches as a list. '''
    
    # Is a branch?
    if type(br) != branch: return list()
    
    # Make list.
    li = list()
    
    # Traverse.
    s = br.child
    li.extend(s.children)
    for bra in s.children: li.extend(getAllBranches(bra))
    return li
        
def isSibling(br,other):
    
    ''' Given a branch, determine if that branch is 
    adjacent to another branch. '''

    if type(br)     != branch: return False
    if type(other)  != branch: return False
    return (br.parent == other.parent)

# Parsing Assets

class parser:
    
    ''' Parsing object for Newick strings representing 
    a phylogenetic tree. '''
    
    def __init__(self, newick):
        self.newick           = newick
        self.parsed_structure = None
        
    def parse(self):
        
        ''' Parse the stored newick string into a topological structure. '''
        
        # Get the newick string.
        read = self.newick
        
        # Check if proper newick format.
        if read[-1] != ';':
            raise ParsingError('Input newick string does not end with ";".')
        
        length   = len(read)
        toplevel = parseNewick(read,0,length-1,None)
        assignParents(toplevel)
        self.parsed_structure = toplevel
        return toplevel
    
    def __str__(self):
        
        if not self.parsed_structure: self.parse()
        return str(self.parsed_structure) + ';'

def getBalancingBracket(newick,i):
    
    ''' Given a position of an opening bracket in a 
    newick string, i, output the closing bracket's
    position that corresponds to this opening bracket. '''
    
    opencount = 0
    while (i+1 < len(newick)):
        if newick[i+1] == '(': opencount += 1 
        elif newick[i+1] == ')':
            opencount -= 1
            if (opencount < 0):
                return i+1
        i += 1
    raise ParsingError('Could not find balancing bracket.')

def getBranchLength(newick,i):
    
    ''' Given a position of a colon symbol (indicating a 
    branch length), return the branch length. '''
    
    runningstr = ''
    while (i+1 < len(newick)):
        if not newick[i+1].isdigit() and not newick[i+1] == '.':
            break
        else: runningstr += newick[i+1]
        i += 1
    try: out = float(runningstr)
    except:
        raise ParsingError('Could not parse branch length as float at %d.' % (i))
    return (i,out)

def getLeafName(newick,i):
    
    ''' Given the position of a leaf, find its complete name. '''
    
    runningstr = ''
    while (i < len(newick)):
        if newick[i] in [',',':',';',')'] \
           or i == len(newick):
            return (i-1,runningstr)
        runningstr += newick[i]
        i += 1
    return (i-1,runningstr)
       
def parseNewick(newick,i,j,top):
    
    ''' Parse a newick string into a topological newick
    structure given a top-level node. '''
    
    # Input is substring [i,j] and top-level node.
    first = None
    while (i <= j):
        
        # Create new node.
        newnode = node()
        
        # Construct branch edge.
        if top != None:
            li = top.children
            nb = branch(newnode,0.0,top)
            li.append(nb)
        elif not first: first = newnode
        
        # Hit a subtree?
        if newick[i] == '(':
            k, _ = (getBalancingBracket(newick,i)), None
            brlength = 0.0
            if newick[k+1] == ':': # Has branch length?
                _, brlength = getBranchLength(newick,k+1)
            elif newick[k+1] != ',': # Has name or statistical support?
                _, name       = getLeafName(newick,k+1)
                newnode.label = name
                if newick[_+1] == ':': # Branch length after this?
                    _, brlength = getBranchLength(newick,_+1)
            if brlength > 0.0 and top: nb.branch_length = brlength
            parseNewick(newick,i+1,k,newnode)
            if (_ != None): k = _
        
        # Hit a leaf?
        else:
            k, name = getLeafName(newick,i)
            brlength = 0.0
            if (k+1 <= j):
                if newick[k+1] == ':':
                    k, brlength = getBranchLength(newick,k+1)
            if brlength > 0.0 and top: nb.branch_length = brlength
            newnode.label = name
            
        # Advance forward.
        i = k + 1 # Next symbol.
        if (i<j):
            while (i+1 <= j and newick[i] != ','): i += 1
        i += 1

    return first