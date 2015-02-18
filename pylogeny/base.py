''' Definitions for generalized containers and objects used by other structures
in this framework. '''

# Date:   Nov 26 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

# Imports

from collections import Sized, Iterable, Container

# Function Definitions

def LCS(str1,str2): return longest_common_substring(str1,str2)

def longest_common_substring(s1,s2):
    
    ''' Simplified, traditional LCS algorithm implementation. '''
    
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else: m[x][y] = 0
    return s1[x_longest - longest: x_longest]

# Tree Elements

class treeStructure(Container):
    
    ''' Defines a base collection of treeNodes and treeBranches in a
    hierarchical tree structure. '''
    
    root = None
    
    def __init__(self,root=None):
        
        if root == None:
            self.root = treeNode()
        else: self.root = root
    
    def __contains__(self,x):
        
        ''' Determines whether a node is found in the tree
        structure. '''
        
        for node in self.getPostOrderTraversal():
            if node == x: return True
        return False
    
    def getRoot(self):
        
        ''' Return the top-level, root, node of 
        the tree. '''
        
        return self.root    
    
    @staticmethod
    def leaves(root):
        if root.isLeaf(): return [root]
        li = list()
        for branch in root.getChildren():
            li.extend(treeStructure.leaves(branch.getChild()))
        return li        
        
    def getAllLeaves(self): 
        
        return treeStructure.leaves(self.root)
    
    @staticmethod
    def nodes(root):
        if root == None: return []
        li = list()
        li.append(root)
        for branch in root.getChildren():
            li.extend(treeStructure.nodes(branch.getChild()))
        return li    
    
    def getAllNodes(self):
        
        return treeStructure.nodes(self.root)
    
    @staticmethod
    def postOrderTraversal(root):
        if (root == None): return []
        li = list()
        for branch in root.getChildren():
            li.extend(treeStructure.postOrderTraversal(branch.child))
        li.append(root)
        return li
        
    def getPostOrderTraversal(self):
        
        return treeStructure.postOrderTraversal(self.root)

    def _createStringLinesForNode(self,root,levelsDeep=0):
        if (root == None): return '\n'
        elif (len(root.getChildren()) == 0):
            return '  '*levelsDeep + '| ' + root.label
        strOut = ''
        for branch in root.getChildren():
            strOut += self._createStringLinesForNode(branch.child,levelsDeep+1)
        return strOut
    
    def __str__(self):
        
        ''' Returns a string representation of the tree. '''

        return self._createStringLinesForNode(self.root)

class treeNode(object):
    
    ''' A node in a tree. '''

    label    = None
    parent   = None
    children = None
    
    def __init__(self,lbl=None,children=None,parent=None):
        
        self.label  = lbl
        self.parent = parent
        if children != None:
            self.children = children
        else: self.children = list()
            
    def getLabel(self): return self.label
    def getParent(self): return self.parent
    def addChild(self,item):
        self.children.append(item)
    def getChildByIndex(self,i): return self.children[i]
    def getChildren(self): return self.children
    def isLeaf(self): return len(self.children) == 0
    def isInternalNode(self): return (len(self.children) > 0)

class treeBranch(object):
    
    ''' A branch in a tree. '''
    
    label  = ''    
    parent = None
    child  = None
    
    def __init__(self,parent=None,child=None,label=''):
        
        self.parent = parent
        self.child  = child
        self.label  = label
        
    def getLabel(self): return self.label
    def getParent(self): return self.parent
    def getChild(self): return self.child  

# Tries

class trieNode(treeNode):
    
    ''' A subclass of treeNode that allows for checking 
    non-zero members amongst children branches and other
    conveniences. '''
    
    def getParentNode(self): return self.parent.parent
    
    def setChildNode(self,child,newchild):
        
        for branch in self.children:
            if branch.child != None and branch.child == child:
                branch.child = newchild
                if newchild: newchild.parent = branch
                return
    
    def iterNonEmptyChildrenNodes(self):    
        
        ''' Iterate over all children that are not empty. '''
        
        for branch in self.children:
            if (branch.child != None):
                yield branch.child
                
    def getNonEmptyChildrenNodes(self):
        
        ''' Acquire a list of all non-empty children. '''
        
        return [x for x in self.iterNonEmptyChildrenNodes()]
    
    def getNonEmptyChildrenBranches(self):
        
        ''' Acquire a list of all non-empty children. '''
    
        return [x.parent for x in self.iterNonEmptyChildrenNodes()]
    
    def getNonEmptyChildrenBranchLabels(self):
        
        return [x.parent.label for x in self.iterNonEmptyChildrenNodes()]
    
    def numEmptyChildrenNodes(self):
        
        ''' Acquire the number of children nodes that are marked 0 
        or nonexistent. '''
        
        empty = 0
        
        for branch in self.children:
            if (branch.child == None or branch.child.label == 0):
                empty += 1
                
        return empty

class trie(Sized,treeStructure):
    
    ''' Defines a trie across a range of strings. '''
    
    alphabet  = None
    root      = None
    count     = 0
    nextLabel = 1
    
    def __init__(self):
        
        self.alphabet = []
        self.root = trieNode(lbl=0)
        
    def __contains__(self,x):

        ''' Implementing for interface (Container). '''
        
        return self.search(x)
        
    def __len__(self):
        
        ''' Implementing for interface (Sized). '''
        
        return self.count
        
    def getAlphabet(self): return self.alphabet
    def getRoot(self): return self.root
    
    def _makeNewChildForNodes(self,node,char):
        
        ''' Called when a new item is added to the alphabet. Adds a new
        empty node to all children. '''
        
        if node == None: return
        if len(node.getChildren()) != 0:
            for child in node.getChildren():
                self._makeNewChildForNodes(child.child,char)
        node.addChild(treeBranch(parent=node,label=char))
    
    def _newNode(self,parent):
        
        ''' Create a new node. '''
        
        t = trieNode(parent=parent,lbl=0)
        for char in self.alphabet:
            t.addChild(treeBranch(parent=t,label=char))
        return t
    
    def _deleteNode(self,node,child=None):
        
        ''' Delete a node. '''
        
        # If a child was given, clear it.
        if child != None:
            node.setChildNode(child,None)
            child.parent = None
            
        # If this is root, stop.
        if node == self.root: return
        
        # Check node for emptiness.
        numempty = node.numEmptyChildrenNodes()
        if (numempty == len(self.alphabet)):
            # All empty.
            self._deleteNode(node.getParentNode(),node)
    
    def search(self,seq):
        
        ''' Search for a sequence in the trie. Returns true
        if it exists. '''
        
        index  = 0
        cursor = self.getRoot()
        
        for character in seq:
            
            # In alphabet?
            if not character in self.alphabet:
                return False
            
            # Get index in alphabet.
            index = self.alphabet.index(character)
            
            # Traverse trie.
            nextitem = cursor.getChildByIndex(index)
            if (nextitem.child == None):
                return False
            cursor = nextitem.child
        
        return (cursor.label > 0)
    
    def insert(self,seq):
        
        ''' Dynamically insert a sequence into the trie. '''
    
        index  = 0
        cursor = self.getRoot()
        self.count += 1
        
        for character in seq:
            
            # In alphabet?
            if not character in self.alphabet:
                self.alphabet.append(character)
                self._makeNewChildForNodes(self.getRoot(),character)
            
            # Get index in alphabet.
            index = self.alphabet.index(character)
            
            # Traverse trie.
            nextitem = cursor.getChildByIndex(index)
            if (nextitem.child == None):
                nextitem.child = self._newNode(nextitem)
            cursor = nextitem.child
        
        # Mark the leaf.
        cursor.label = self.nextLabel
        self.nextLabel += 1
        
        # Return its label.
        return cursor.label
        
    def delete(self,seq):
        
        ''' Remove a sequence from the trie. Will not remove
        added characters to alphabet. '''
        
        index  = 0
        cursor = self.getRoot()
        
        for character in seq:
            
            if not character in self.alphabet:
                raise LookupError('Sequence %s not in trie.' % (seq))
            index = self.alphabet.index(character)
            
            # Traverse trie.
            nextitem = cursor.getChildByIndex(index)
            if (nextitem.child == None):
                raise LookupError('Sequence %s not in trie.' % (seq))
            cursor = nextitem.child
        
        # Work way up and check for completely empty subtrees.
        self.count -= 1
        self._deleteNode(cursor)
        
    def _createStringLinesForNode(self,root,levelsDeep=-1):
        
        ''' Private method for stringifying this object. '''
        
        if (root == None): return ''
        elif (root != self.root):
            strOut   = '  '*levelsDeep + '| ' + root.parent.label + '\n'
        else: strOut = ''
        for branch in root.getChildren():
            strOut += self._createStringLinesForNode(branch.child,levelsDeep+1)
        return strOut     

class patriciaTree(trie):
    
    ''' Defines a PATRICIA tree (condensed trie) across a range of strings. '''
    
    def __contains__(self,x): return (self.search(x) != 0)
    
    def _deleteNode(self,node,child=None):
        
        ''' Merges node upward if it only has 1 child. '''
        
        # If a child was given, clear it.
        if child is not None:
            node.setChildNode(child,None)
            child.parent = None
            
        # If this is root, stop.
        if node is self.root: return
        
        # Check node for emptiness.
        numEmpty = node.numEmptyChildrenNodes()
        alphaLen = len(self.alphabet)
        if (numEmpty == alphaLen - 1):
            # Only has a single child. Merge upwards.
            parent     = node.getParentNode()
            edge       = node.parent
            nonemptybr = node.getNonEmptyChildrenBranches()[0]
            edge.label = edge.label + nonemptybr.label
            parent.setChildNode(node,nonemptybr.child)    
        elif (numEmpty == alphaLen):        
            # All empty.
            self._deleteNode(node.getParentNode(),node)
    
    def _query(self,seq):
        
        ''' Acquires the leaf for a sequence if it exists. Returns the
        node and its parent of wherever failure or success occurs. Returns
        a tuple of elements where cursor is the node the query stopped at,
        parentofcursor is the branch associated with that node, cpos is the
        position of seq that the query stopped at, and spos is where cpos was
        last at before finding a mismatch. '''
        
        index,cpos,spos = 0,0,0
        cursor = self.getRoot()
        parentofcursor = cursor
        seqLength = len(seq)
        
        while (cpos < seqLength):

            # Get character.
            character = seq[cpos]
            if not character in self.alphabet:
                return None, None, cpos, spos
            index = self.alphabet.index(character)
            
            # Go to next respective item in trie.
            nextitem       = cursor.getChildByIndex(index)
            parentofcursor = nextitem
            if (nextitem.child == None):
                # Item is not in the trie.
                return None, parentofcursor, cpos, spos
            cursor = nextitem.child 
            label  = nextitem.label
            spos   = cpos
            lpos   = 0
            
            # Check character matching.
            labelLength = len(label)
            while (lpos < labelLength):
                if (cpos >= seqLength):
                    # Still more in the label! Not found.
                    return None, parentofcursor, cpos, spos
                elif (label[lpos] != seq[cpos]):
                    # Found a mismatch.
                    return None, parentofcursor, cpos, spos
                lpos += 1
                cpos += 1
                
        return cursor, parentofcursor, seqLength, spos
    
    def search(self,seq):
        
        ''' Search for a sequence in the PATRICIA tree. Returns its position in
        addition sequence if it exists. Else, returns 0. '''
        
        query,_,_,_ = self._query(seq)
        if (query is None): return 0
        return query.label       
    
    def insert(self,seq):
        
        ''' Dynamically insert a sequence into the PATRICIA tree. Returns the
        unique index in the tree for that string. '''
        
        # Check alphabet membership.
        for character in seq:
            if not character in self.alphabet:
                self.alphabet.append(character)
                self._makeNewChildForNodes(self.getRoot(),character)        
        
        # Perform a query.
        query,queryparent,cpos,spos = self._query(seq)
        if query == None or query.label == 0: self.count += 1
        else: return 0 # Already in the tree.
        
        prevLabel = ''
        # Insert Case 1: Middle of Edge. Check for label of edge.
        currEdge = queryparent
        cmpLabel = seq[spos:cpos+1]
        
        if (len(currEdge.label) > 1 and cmpLabel != currEdge.label):
            
            # Know point of convergence.
            newLabel = seq[spos:cpos]
            oldLabel = currEdge.label
            
            # Replace previous edge label.
            currEdge.label = newLabel
            
            # See what has to be split to a new edge.
            extent = oldLabel[len(newLabel):]
            
            # Make new node to join split edge.
            oldNode = currEdge.child
            newNode = self._newNode(currEdge)
            currEdge.child = newNode
            
            # Get ready for insertion to newChild.
            extentIndex  = self.alphabet.index(extent[0])
            extentBranch = newNode.getChildByIndex(extentIndex)            
            
            # Special Case: Adding a subset of an existing entry.
            if cpos >= len(seq):
                # Simply mark the new node.
                newNode.label = self.nextLabel
                insertBranch  = None
            
            # Add extent to newChild as well as this insertion.
            else:
                insertIndex  = self.alphabet.index(seq[cpos])
                insertBranch = newNode.getChildByIndex(insertIndex)
            extentBranch.label = extent
            extentBranch.child = oldNode
            oldNode.parent = extentBranch            
            
        # Insert Case 2: New Node.
        else: insertBranch = queryparent
        
        if insertBranch != None:
            prevLabel = insertBranch.label
            if cpos != len(seq):
                insertBranch.label = seq[cpos:]
            if insertBranch.child is None:
                insertBranch.child = self._newNode(insertBranch)
            insertBranch.child.label = self.nextLabel
        
        self.nextLabel += 1
        return self.nextLabel-1

    def delete(self,seq):
        
        ''' Remove a sequence from the PATRICIA tree. Will not remove 
        added characters to alphabet. '''

        query,_,_,_ = self._query(seq)
        
        if query is None: raise LookupError('Sequence %s not in trie.' % (seq))
        
        # Work way up and check for completely empty subtrees.
        self.count -= 1
        self._deleteNode(query)
