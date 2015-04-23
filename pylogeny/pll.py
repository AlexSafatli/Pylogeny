''' Wrap C extension for libpll library for use in natural Python. '''

# Date:   Feb 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import alignment, rearrangement
from libpllWrapper import *
from tempfile import NamedTemporaryFile as NTempFile
import os, sys

class dataModel:
    
    ''' Encapsulating a phylogenetic tree (as topology) + corresponding alignment
    into a libpll-associated data structure. Allows for log-likelihood
    scoring of this model. MUST BE CLOSED AFTER USE. '''

    def __init__(self,topo,alignm,model=None):

        ''' Initialize the data model and respective structures. 
        
        :param topo: a topology object
        :type topo: :class:`.rearrangement.topology`
        :param alignm: a phylip-friendly alignment object.
        :type alignm: :class:`.alignment.phylipFriendlyAlignment`
        
        '''
        
        # Ensure that given alignment object is of correct type.
        if type(alignm) != alignment.phylipFriendlyAlignment:
            raise TypeError('libpll interfacing needs phylipFriendlyAlignment.')

        # Model handling. See if one already defined in alignment.
        if not model:
            if hasattr(alignm,'_pllmodel'): self.model = alignm._pllmodel
            else:
                self.model = partitionModel(alignm)
                self.model.createSimpleModel()
                alignm._pllmodel = self.model
                alignm.paths['pll'] = self.model.getFileName()
        else: self.model = model        

        # Other things.
        modf           = self.model.getFileName()
        self.instance  = new(alignm.getPhylip(),
                             topo.toUnrootedNewick(),
                             modf)
        
    def getNewickString(self):
        
        ''' Acquire the Newick string of the problem instance.
        
        :return: a Newick string
        
        '''
        
        return getNewickString(self.instance)
        
    def getLogLikelihood(self):
        
        ''' Calculates log-likelihood using libpll.
        
        :return: a floating point value
        
        '''
       
        return getLogLikelihood(self.instance)

    def close(self):

        ''' If done with this particular problem. Frees associated memory. '''
        
        destroy(self.instance)

class partitionModel:
    
    ''' A partition model intended for libpll. '''
    
    def __init__(self,ali):
        
        ''' Initialize a partition model (for internal use by libpll). 

        :param ali: an alignment object
        :type ali: :class:`.alignment.alignment`

        '''
        
        self.handle = NTempFile(delete=False)
        self.length = len(ali)
        self.isprotein = (ali.getDataType() == 'protein')
        
    def getFileName(self):
        
        ''' Get the file name of the model file.
        
        :return: a string
        
        '''
        
        return self.handle.name
    
    def createSimpleModel(self,pmodel='WAG'):
        
        ''' Establish a simple model (e.g., one type).
        
        :param pmodel: optional; what protein model to use (as described in pll)
        :type pmodel: a string (default 'WAG')
        :return: None
        
        '''
        
        if self.isprotein: simplemodel = "%s, p1 = 1-%d\n" % (
            pmodel,self.length)
        else: simplemodel = "DNA, p1 = 1-%d\n" % (self.length)
        self.handle.write(simplemodel)
        self.handle.close()
        
    def createModel(self,models,partnames,ranges):
        
        ''' Establish a more complex model.
        
        :param models: a list of model names (e.g., 'WAG', 'DNA')
        :type models: a list of strings
        :param partnames: a list of partition names (e.g., 'p1', 'p2')
        :type partnames: a list of strings
        :param ranges: a list of range tuples (what ranges of alignment)
        :type ranges: a list of integer tuples
        :return: None
        
        '''
        
        for m in xrange(len(models)):
            mod    = models[m]
            par    = partnames[m]
            rl, ru = ranges[m]
            towr   = "%s, %s = %d-%d\n" % (mod,par,rl,ru)
            self.handle.write(towr)
        self.handle.close()
        
    def close(self):
        
        ''' Delete file. '''
        
        os.unlink(self.handle.name)
