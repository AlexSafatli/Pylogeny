''' Wrap C extension for libpll library for use in natural Python. '''

# Date:   Feb 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import alignment, rearrangement
from libpllWrapper import *
from tempfile import NamedTemporaryFile as NTempFile
from cStringIO import StringIO
import os, sys

class dataModel:
    
    ''' Encapsulating a phylogenetic tree (as topology) + corresponding alignment
    into a libpll-associated data structure. Allows for log-likelihood
    scoring of this model. **MUST BE CLOSED AFTER USE.** '''

    def __init__(self,topo,alignm,model=None):

        ''' Initialize the data model and respective structures. 
        
        :param topo: A topology object.
        :type topo: :class: `rearrangement.topology`
        :param alignm: A phylipFriendlyAlignment object.
        :type alignm: :class: `alignment.phylipFriendlyAlignment`
        
        '''
        
        # Ensure that given alignment object is of correct type.
        if type(alignm) != alignment.phylipFriendlyAlignment:
            raise TypeError('libpll interfacing needs phylipFriendlyAlignment.')
        
        # Declarations
        isProtein = (alignm.data.dataType == 'protein')

        # Model handling. See if one already defined in alignment.
        if not model:
            if hasattr(alignm,'_pllmodel'): self.model = alignm._pllmodel
            else:
                self.model = partitionModel(alignm)
                self.model.createSimpleModel(isProtein)
                alignm._pllmodel = self.model
                alignm.paths['pll'] = self.model.getFileName()
        else: self.model = model        

        # Other things.
        modf           = self.model.getFileName()
        self.instance  = new(alignm.getPhylip(),
                             topo.toUnrootedNewick(),
                             modf)
        
    def getNewickString(self):
        
        ''' Acquire the Newick string of the problem instance. '''
        
        return getNewickString(self.instance)
        
    def getLogLikelihood(self):
        
        ''' Calculates log-likelihood using libpll. '''
       
        return getLogLikelihood(self.instance)

    def close(self):

        ''' If done with this particular problem. Frees associated memory. '''
        
        destroy(self.instance)

class partitionModel:
    
    ''' A partition model intended for libpll. '''
    
    def __init__(self,ali):
        
        self.handle = NTempFile(delete=False)
        self.length = len(ali)
        
    def getFileName(self):
        
        ''' Get the file name of the model file. '''
        return self.handle.name
    
    def createSimpleModel(self,protein,pmodel='WAG'):
        
        ''' Establish a simple model (e.g., one type). '''
        
        if protein: simplemodel = "%s, p1 = 1-%d\n" % (pmodel,self.length)
        else:       simplemodel = "DNA, p1 = 1-%d\n" % (self.length)
        self.handle.write(simplemodel)
        self.handle.close()
        
    def createModel(self,models,partnames,ranges):
        
        ''' Establish a more complex model. '''
        
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
