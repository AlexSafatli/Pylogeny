''' Toolkit for performance of Parsimonious Criterion (Parsimony) methods of
optimization of a phylogenetic topology with a particular set of data. '''

# Date:   Mar 3 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from base import treeStructure
postorder = treeStructure.postOrderTraversal

class profile_set:
    
    ''' Hold a set of site_profile profiles for an 
    entire alignment. '''

    def __init__(self,alignment):
        
        self.alignment = alignment
        self.numSites  = alignment.getSize()
        self.taxa      = {}
        self.profiles  = []
        self.weights   = []
        self.sites     = []
        self._constructSet()
        self._buildTaxaDict()
    
    def _constructSet(self):
        
        ''' PRIVATE: Build all site_profiles. '''
        
        num = self.numSites
        for site in xrange(num):
            pro = site_profile(self.alignment,site)
            self.sites.append(pro)
            if pro in self.profiles:
                self.weights[self.profiles.index(pro)] += 1
            else:
                self.profiles.append(pro)
                self.weights.append(1)
                
    def _buildTaxaDict(self):
        
        ''' PRIVATE: Build taxa dictionary. '''
        
        tax = self.alignment.getTaxa()
        for t in xrange(len(tax)): self.taxa[tax[t]] = t
                
    def __len__(self):
        
        return len(self.profiles)
    
    def weight(self,val):
        
        return self.weights[val]
    
    def get(self,val):
        
        return self.profiles[val]
    
    def getForTaxa(self,val,tax):
        
        return self.profiles[val].vector[self.taxa[tax]]

class site_profile:
    
    ''' Consolidate the single-column alignment at
    a region into a set of components on the basis
    of similarity alone. '''
    
    def __init__(self,alignment,site):
        
        self.alignment = alignment
        self.site      = site
        self.alphabet  = None
        self.vector    = ''
        self._buildVector()
        
    def __eq__(self,o):
        
        if (o == None): return False
        return (self.vector == o.vector)
    
    def __ne__(self,o):
        
        return not self.__eq__(o)
        
    def __str__(self):
        
        return self.vector
        
    def _buildVector(self):
        
        ''' Build the vector comprising the set of components. '''
        
        aliData = self.alignment.data
        sitesli = aliData.sequenceSlice(self.site)
        compn   = 0
        comps   = {}
        for it in sitesli:
            if not it in comps:
                comps[it] = compn
                compn += 1
            self.vector += str(comps[it])
        self.alphabet = comps

def fitch_cost(topology,profiles):
    
    ''' Calculate the cost using Fitch algorithm on 
    profile set and alignment. Deprecated: Python implementation 
    of the Fitch algorithm; see fitch C++ module for 
    a C++ implementation that is roughly four times faster. '''
    
    # Calculate parsimony score.
    total = 0
    porde = [x for x in postorder(topology.getRoot())]
    for profile in xrange(len(profiles)):
        posdata, local = {}, 0
        for node in porde:
            # If is leaf.
            if (len(node.children) == 0):
                posdata[node] = profiles.getForTaxa(profile,node.label)
            # Not a leaf.
            else:
                a = posdata[node.children[0].child]
                b = posdata[node.children[1].child]
                # Get intersection.
                X = set(a).intersection(set(b))
                # If intersection is nothing...
                if (len(X) == 0):
                    # Cost goes up, take union.
                    local += 1
                    X = set(a).union(set(b))
                posdata[node] = X
        total += local*profiles.weight(profile)
    return total

def fitch(topology,alignment):
        
    ''' Perform the Fitch algorithm on a given
    tree topology and associated alignment. Deprecated: 
    Python implementation of the Fitch algorithm; 
    see fitch C++ module for a C++ implementation 
    that is roughly four times faster. '''
    
    # Turn alignment sets into profiles.
    profiles = profile_set(alignment)
    if (len(profiles) == 0): return None
    return fitch_cost(topology,profiles)