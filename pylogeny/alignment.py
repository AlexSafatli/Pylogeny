''' Handle input biological sequence alignment files for the purposes of phylogenetic inference. Will read all types of alignment files by utilizing the P4 python phylogenetic library. '''

# Date:   Jan 24 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

# Imports

import os, p4, model, tree, base
from newick import newickParser
from executable import fasttree
from tempfile import NamedTemporaryFile as NTempFile
from shutil import copyfile
p4.var.warnReadNoFile = False

# Class Definitions

class alignment(object):
    
    ''' Wrap a biological sequence alignment to enable functionality
    necessary for phylogenetic inference. Makes use of temporary
    files; requires to be closed once no longer needed. '''
    
    def __init__(self,inal=None):
        
        ''' Instantiate an object intended to wrap an alignment for the purposes of running phylogenetic inference. 
        
        :param inal: An alignment file path (most formats are accepted). 
        
        '''
        
        # Intergrate the p4 phylogenetic library.
        if inal is None:
            self.data = p4.Alignment()
            
        else:
            p4.read(inal) # Read the alignment file/string.
            self.data  = p4.var.alignments[-1]

        # Augment alignment with a discrete state model.
        self.model = model.DiscreteStateModel(self)

        # Keep track of temporary files.
        self.paths = {} 

    # Internals

    def __getitem__(self,i): return self.data.sequences[i]
    def __str__(self):       return str(self.data)
    def __len__(self):       return len(self.data)
    
    def __iter__(self):
        
        for i in xrange(self.getNumSeqs()): yield self[i]

    def _makeFASTAFile(self):
        
        ''' PRIVATE: Make a temporary FASTA file. '''
        
        if ('fasta' in self.paths): return
        fh = NTempFile(delete=False)
        self.paths['fasta'] = fh.name
        self.data.writeFasta(fh.name)
        fh.close()

    def close(self):
        
        ''' Delete all temporary files and clear data. '''
        
        self.data  = None
        self.model = None
        if len(self.paths) > 0:
            for f in self.paths.values(): os.unlink(f)
            self.paths = {}

    # Accessors
    
    def toStrList(self):
        
        ''' Get all sequences as a list of strings. '''
        
        return [self.getSequence(i) for i in xrange(self.getNumSeqs())]
    
    def getStateModel(self): return self.model
    
    def getSize(self):
        
        ''' Return the size of the alignment, or how many characters
        there are in each respective item in the alignment. '''
        
        return self.__len__()
    
    def getNumSeqs(self):
        
        ''' Return the number of sequences that are present in the 
        sequence alignment. '''
        
        return len(self.data.sequences)
    
    def getDim(self):
        
        ''' Return the dimensionality of the sequence alignment (how 
        many different types of characters). '''
        
        return self.data.dim
    
    def getSequence(self,i):
        
        ''' Acquire the ith sequence. '''
        
        return self[i].sequence

    def getFASTA(self):
        
        ''' Get (and create if not already) a path to a temporary 
        FASTA file. This will be deleted upon closure of the alignment
        instance. '''
        
        self._makeFASTAFile()
        return self.paths['fasta']

    def getApproxMLNewick(self):
        
        ''' Get a tree in newick format via use of FastTree that serves as
        an approximation of the maximum likelihood tree for this data. '''
        
        ft = fasttree(
            self.getFASTA(),isProtein=(self.data.dataType=='protein'))
        return ft.run()
    
    def getApproxMLTree(self):
        
        ''' Get a tree object for an approximation of the maximum likelihood
        tree for this data using FastTree. '''
        
        return tree.tree(self.getApproxMLNewick(),check=True)

    def getFastTreeNewick(self):
        
        ''' Alias for the "getApproxMLNewick()" function. '''
        
        return self.getApproxMLNewick()

    def getTaxa(self):
        
        ''' Return taxa names. '''
        
        return self.data.taxNames

    # Interface Revealing
    
    def getAlignment(self):
        
        ''' Acquire the alignment data structure (P4 module). '''
        
        return self.data
    
    def bootstrap(self):
        
        ''' Perform bootstrapping on the alignment data. '''
        
        self.data = self.data.bootstrap()    


class phylipFriendlyAlignment(alignment):
    
    ''' An alignment object that renames all comprising taxa in order
    to be able to be written as a strict Phylip file. '''
    
    def __init__(self,inal=None):
        
        # Call superclass constructor.
        super(phylipFriendlyAlignment,self).__init__(inal)
        
        # Keep track of all taxa names.
        self.namedict = {}

        # Acquire a temporary file handle for Nexus with proper names.
        self._makeProperNexusFile()

        # Reassign all taxa names to shorter ones (for Phylip).
        self._reassignNames()  

    def _makePhylipFile(self):
        
        ''' PRIVATE: Make a temporary Phylip file. '''
        
        if ('phylip' in self.paths): return
        fh = NTempFile(delete=False)
        self.paths['phylip'] = fh.name
        self.data.writePhylip(
            fh.name,whitespaceSeparatesNames=False)
        fh.close()

    def _makeProperNexusFile(self):
        
        ''' PRIVATE: Make a temporary Nexus file. '''
    
        if ('nexus' in self.paths): return
        fh = NTempFile(delete=False)
        self.paths['nexus'] = fh.name
        self.data.writeNexus(
            fName=fh.name)    
        fh.close()

    def _reassignNames(self):
        
        ''' PRIVATE: Reassign all sequence names to shorter 
        Phylip-friendly integer-based names. '''
        
        d  = self.namedict
        for item in self.data.sequences:
            _         = item.name
            ind       = self.data.sequences.index(item)
            name      = 'T' + str(ind)
            d[name]   = _
            item.name = name
            self.data.taxNames[ind] = name

    def getPhylip(self):
        
        ''' Get a path to a temporary Phylip file. This will 
        be deleted upon closure of the alignment instance. '''
        
        self._makePhylipFile()
        return self.paths['phylip']

    def writeProperNexus(self,wri):
        
        ''' Write a Nexus file with proper names. '''
        
        copyfile(self.paths['nexus'],wri)
        return True
    
    def reassignFromReinterpretedNewick(self,tr):
        
        ''' Replace all proper names with reassigned names in a Newick tree. '''
        
        p = newickParser(tr).parse()
        nodes = base.treeStructure(p).getAllNodes()
        for node in nodes:
            l = node.label
            if l in self.namedict.values():
                for n in self.namedict:
                    if self.namedict[n] == l:
                        node.label = n
                        break
        return str(p) + ';'
    
    def reinterpretNewick(self,tr):
        
        ''' Replaces all reassigned names to proper names in a Newick tree. '''
               
        p = newickParser(tr).parse()
        nodes = base.treeStructure(p).getAllNodes()
        for node in nodes:
            l = node.label
            if l in self.namedict: node.label = self.namedict[l]
        return str(p) + ';'
    
    def getProperName(self,n):
        
        ''' Return the actual name for an integer-based sequence name 
        that was reassigned at initialization. '''
        
        if n in self.namedict: return self.namedict[n]
        else: return None
    
    def getTaxa(self):
        
        ''' Return current taxa names in the alignment. '''
        
        return [x for x in self.namedict]
    
    def recreateObject(self):
        
        ''' Reintializes the object. '''
        
        return phylipFriendlyAlignment(self.getFASTA())