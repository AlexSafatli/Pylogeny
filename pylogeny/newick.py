''' Newick string parsing and object interaction. A Newick string can represent
a phylogenetic tree. '''

# Date:   Feb 16 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from rearrangement import topology
from random import shuffle
from math import factorial as fact
import base

# Exception Handling

class ParsingError(Exception):
    def __init__(self,val): self.value = val
    def __str__(self): return repr(self.value)

# Class Definitions

class node(base.treeNode):
    
    ''' Node for a tree parsed from a Newick string. '''
        
    def __init__(self,lbl='',children=None,parent=None):
        
        ''' Initialize a node in a tree parsed from a Newick string.
        
        :param lbl: a label for this node
        :type lbl: a string
        :param children: an optional set of children (branches); default none
        :type children: a list of :class:`.branch` objects
        :param parent: an optional parent branch for this node; default none
        :type parent: a :class:`.branch` object
        
        '''
        
        super(node,self).__init__(lbl,children,parent)
        
    def __str__(self):
        sT, sl = self.children, self.label
        if len(self.children) > 0:
            # Perform sorting to ensure consistent naming scheme.
            l = sorted([st for st in sT],key=lambda d: d.child.label)
            return '(%s)%s' % (','.join(map(str,l)),sl)
        else: return self.label

class branch(base.treeBranch):
    
    ''' Branch for a tree parsed from a Newick string. '''
    
    def __init__(self,chi,l,parent=None,s=None):
        
        ''' Initialize a branch in a tree parsed from a Newick string. 
        
        :param chi: a child node
        :type chi: a :class:`.node` object
        :param l: a branch length
        :type l: a floating point value
        :param parent: an optional parent node; default none
        :type parent: a :class:`.node` object
        
        '''
        
        self.parent         = parent
        self.child          = chi
        self.branch_length  = l
        self.branch_support = s
        
    def __str__(self):
        st, sl = self.child, self.branch_length
        if sl > 0: return '%s:%s' % (str(st),str(sl)) # Branch Length
        else:      return '%s'    % (str(st))         # No Branch Length

# Traversal Functions

def assignParents(top):
    
    ''' Should be a one-time use function. Goes through
    and assigns parents to the parsed newick tree structure
    nodes and branches to allow for up-traversal.
    
    :param top: a top-level node for a tree (root node)
    :type top: a :class:`.node` object
    
    '''
    
    # Is a node?
    if type(top) != node: return False
    
    # Traverse.
    for br in top.children:
        br.parent       = top
        br.child.parent = br
        assignParents(br.child)
    return True

def removeBranchLengths(top):
    
    ''' Goes through and removes any stored branch lengths.
    
    :param top: a top-level node for a tree (root node)
    :type top: a :class:`.node` object
    
    '''

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
    node. 
    
    :param top: a top-level node for a tree (root node)
    :type top: a :class:`.node` object    
        
    '''
    
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
    function for rerooting a tree.
    
    :param target: a target node
    :type target: a :class:`.node` object
    :param top: a top-level node for a tree (root node)
    :type top: a :class:`.node` object
    
    '''
        
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

def shuffleLeaves(top):
    
    ''' DANGEROUS: Given a top-level node, shuffle all 
    leaves in this tree.
    
    :param top: a top-level node for a tree (root node)
    :type top: a :class:`.node` object
    
    '''

    leaves = base.treeStructure.leaves(top)
    names  = [x.label for x in leaves]
    shuffle(names)
    for l in xrange(len(leaves)):
        leaves[l].label = names[l]

def getAllBranches(br):
    
    ''' Given a branch, traverse subtree and return 
    comprising branches as a list.
    
    :param br: a branch from a tree
    :type br: a :class:`.branch` object
    
    '''
    
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
    adjacent to another branch. 
    
    :param br: a branch from a tree
    :type br: a :class:`.branch` object
    :param other: another branch from a tree
    :type other: a :class:`.branch` object
    
    '''

    if type(br)     != branch: return False
    if type(other)  != branch: return False
    return (br.parent == other.parent)

# Parsing Assets

class newickParser:
    
    ''' Parsing object for Newick strings. '''
    
    def __init__(self, newick):
        
        ''' Initialize this parser (with a Newick string). 
        
        :param newick: a Newick string
        :type newick: a string
        
        '''
                
        if newick == None or type(newick) != str:
            raise ValueError('Input <%s> is not a valid string!' %
                             (repr(newick)))
        self.newick           = newick
        self.parsed_structure = None
        
    def parse(self):
        
        ''' Parse the stored newick string into a topological structure.
        
        :return: the top-level root :class:`.node` object
        
        '''
        
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
    position that corresponds to this opening bracket.
    
    :param newick: a Newick string
    :type newick: a string
    :param i: a position in the string (index)
    :type i: an integer < length of the string
    :return: an integer
    
    '''
    
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
    branch length), return the branch length.
    
    :param newick: a Newick string
    :type newick: a string
    :param i: a position in the string (index)
    :type i: an integer < length of the string
    :return: an integer
    
    '''
    
    runningstr = ''
    while (i+1 < len(newick)):
        if not newick[i+1].isdigit() and not newick[i+1] == '.':
            break
        else: runningstr += newick[i+1]
        i += 1
    try: out = float(runningstr)
    except:
        raise ParsingError('Could not parse branch length as float at %d in %s.'
                           % (i,newick))
    return (i,out)

def getLeafName(newick,i):
    
    ''' Given the position of a leaf, find its complete name.
    
    :param newick: a Newick string
    :type newick: a string
    :param i: a position in the string (index)
    :type i: an integer < length of the string
    :return: an integer
    
    '''
    
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
    structure given a top-level node.
    
    :param newick: a Newick string
    :type newick: a string
    :param i: a starting position to start parsing
    :type i: an integer
    :param j: an end position to stop parsing
    :type j: an integer
    :param top: a top-level node; start parsing with None
    :type top: a :class:`.node` object
    
    '''
    
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