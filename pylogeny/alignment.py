''' Object model defining a sequence alignment (DNA, RNA, protein sequences). 
Handle input biological sequence alignment files for the purposes of 
phylogenetic inference. Will read all types of alignment files by utilizing the 
P4 python phylogenetic library. '''

# Date:   Jan 24 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

# Imports

import os
import p4
import model
import tree
import base
from newick import newickParser
from executable import fasttree
from tempfile import NamedTemporaryFile as NTempFile
from shutil import copyfile

# Set up.

p4.var.warnReadNoFile = False

# Class Definitions


class alignment(object):
    
    ''' Wrap a biological sequence alignment to enable functionality
    necessary for phylogenetic inference. Makes use of temporary
    files; requires to be closed once no longer needed. '''
    
    def __init__(self,inal=None):
        
        ''' Instantiate an object intended to wrap an alignment for the
        purposes of running phylogenetic inference. 
        
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

    def __getitem__(self,i):
        
        return self.data.sequences[i]
    
    def __str__(self):
        
        return '\n'.join(self.toStrList())
    
    def __len__(self):
        
        return len(self.data)
    
    def __iter__(self):
        
        for i in xrange(self.getNumSeqs()): yield self[i]

    def _makeFASTAFile(self):
        
        ''' PRIVATE: Make a temporary FASTA file. '''
        
        if ('fasta' in self.paths):
            return
        fh = NTempFile(delete=False)
        self.paths['fasta'] = fh.name
        self.data.writeFasta(fh.name)
        fh.close()

    def close(self):
        
        ''' Forcefully delete all temporary files and clear data. '''
        
        self.data  = None
        self.model = None
        if len(self.paths) > 0:
            for f in self.paths.values(): os.unlink(f)
            self.paths = {}

    # Accessors

    def toStrList(self):

        ''' Get all sequences as a list of strings. 

        :return: a list of strings
        
        '''

        return [self.getSequenceString(i) for i in xrange(self.getNumSeqs())]

    def getDataType(self):
        
        ''' Get the data type associated with this alignment (e.g., protein). 
        
        :return: a string indicating the data type ('protein', 'DNA')
        
        '''

        return self.data.dataType

    def getStateModel(self):

        ''' Get the state model associated with this alignment. See model
        module for more information.
        
        :return: a :class:`.model.DiscreteStateModel` object
        
        '''

        return self.model

    def getSize(self):

        ''' Return the size of the alignment, or how many characters
        there are in each respective item in the alignment. 
        
        :return: an integer
        
        '''

        return len(self)

    def getNumSeqs(self):

        ''' Return the number of sequences that are present in the 
        sequence alignment. 
        
        :return: an integer
        
        '''

        return len(self.data.sequences)

    def getDim(self):

        ''' Return the dimensionality of the sequence alignment (how 
        many different types of characters). 
        
        :return: an integer
        
        '''

        return self.data.dim

    def getSequenceString(self,i):

        ''' Acquire the ith sequence as a string. 
        
        :param i: an index in the alignment (associated with a sequence)
        :type i: an integer
        :return: a string associated with the sequence
        
        '''

        return self[i].sequence

    def getFASTA(self):

        ''' Get (and create if not already) a path to a temporary 
        FASTA file. This will be deleted upon closure of the alignment
        instance. 
        
        :return: a string associated with a path in the file system
        
        '''

        self._makeFASTAFile()
        return self.paths['fasta']

    def getApproxMLNewick(self):

        ''' Get a tree in newick format via use of FastTree that serves as
        an approximation of the maximum likelihood tree for this data. 
        
        :return: a Newick or New Hampshire string
        
        '''

        ft = fasttree(
            self.getFASTA(), isProtein=(self.data.dataType=='protein'))
        return ft.run()

    def getApproxMLTree(self):

        ''' Get a tree object for an approximation of the maximum likelihood
        tree for this data using FastTree. 
        
        :return: a :class:`.tree.tree` object
        
        '''

        return tree.tree(self.getApproxMLNewick(),check=True)

    def getTaxa(self):

        ''' Get a list of taxa names associated with the alignment. 
        
        :return: a list of strings
        
        '''

        return self.data.taxNames


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

        ''' Get (and create if not already) a path to a temporary 
        Phylip file. This will be deleted upon closure of the alignment
        instance. 

        :return: a string associated with a path in the file system

        '''

        self._makePhylipFile()
        return self.paths['phylip']

    def writeProperNexus(self,wri):

        ''' Write a Nexus file with proper names. 
        
        :param wri: a path to a (existent or unexistent) file to write to
        :type wri: a string
        
        '''

        copyfile(self.paths['nexus'],wri)
        return True

    def convertOriginalNewick(self,tr):

        ''' Return a Newick string with (original) taxa names that are replaced
        with the shortened forms as they are defined in this object.
        
        :param tr: a Newick string
        :type tr: a string
        :return: a Newick string with all replaced names
        
        '''

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

        ''' Revert the replacing of taxa names with shortened names by changing
        them back to their original form.
        
        :param tr: a Newick string
        :type tr: a string
        :return: a Newick string with all replaced names
        
        '''

        p = newickParser(tr).parse()
        nodes = base.treeStructure(p).getAllNodes()
        for node in nodes:
            l = node.label
            if l in self.namedict: node.label = self.namedict[l]
        return str(p) + ';'

    def getProperName(self,n):

        ''' Return the actual name for an integer-based sequence name 
        that was reassigned at initialization. 
        
        :param n: a shortened taxon name from this object
        :type n: a string
        :return: a string (replaced with the original taxon name)
        
        '''

        if n in self.namedict: return self.namedict[n]
        else: return None

    def getTaxa(self):

        ''' Return current taxa names in the alignment.
        
        :return: a list of shortened taxa names
        
        '''

        return [x for x in self.namedict]

    def recreateObject(self):

        ''' Reintializes the object. '''

        return phylipFriendlyAlignment(self.getFASTA())
