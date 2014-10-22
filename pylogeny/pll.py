''' C Extension to wrap libpll library. '''

# Date:   Feb 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from pylibpll import *
from tempfile import NamedTemporaryFile as NTempFile
import os

class dataModel:
    
    ''' Encapsulating a phylogenetic tree (as topology) + corresponding alignment
    into a libpll-associated data structure. Allows for log-likelihood
    scoring of this model. MUST BE CLOSED AFTER USE. '''

    def __init__(self,topo,alignm,model=None):

        ''' Initialize all structures. '''
        
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
        
    def getLogLikelihood(self):
        
        ''' Calculates log-likelihood using libpll. '''
        
        return getLogLikelihood(self.instance)

    def close(self):

        ''' If done with this particular problem. '''
        
        destroy(self.instance)
  

class partitionModel:
    
    ''' A partition model intended for libpll. '''
    
    def __init__(self,ali):
        
        self.handle = NTempFile(delete=False)
        self.length = len(ali)
        
    def getFileName(self):
        
        ''' Get the file name of the model file. '''
        return self.handle.name
    
    def createSimpleModel(self,protein):
        
        ''' Establish a simple model (e.g., one type). '''
        
        if protein: simplemodel = "WAG, p1 = 1-%d\n" % (self.length)
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