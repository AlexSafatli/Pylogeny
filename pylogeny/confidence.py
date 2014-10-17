''' Determine the confidence of a set of trees found in a landscape. Uses TREE-PUZZLE and CONSEL. '''

# Date:   Aug 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import os
from subprocess import PIPE, Popen as system

# Constants
CONSEL_TMPDIR  = '/tmp/consel__py/'
TREEPUZZLE_BIN = 'puzzle'

def AUtest(ls,name):
    c = consel(ls,name)
    c.compute()
    return c.interval

class aTemporaryDirectory(object):
    
    ''' A class intended to be used as a context manager that allows
    Python to run in another directory temporarily. '''
    
    def __init__(self,path=None):
        
        if path == None:
            if not os.path.isdir(CONSEL_TMPDIR):
                os.mkdir(CONSEL_TMPDIR)
            self.path = CONSEL_TMPDIR
        else: self.path = path
        self.current = os.getcwd()
        
    def __enter__(self):      os.chdir(self.path)
    def __exit__(self,*args): os.chdir(self.current)

class treepuzzle(object): # TODO: Move over to executable interface.
    
    ''' Wrap TREE-PUZZLE in order to create an intermediate file for 
    CONSEL to read and assign confidence to a set of trees from a 
    phylogenetic landscape.'''
    
    def __init__(self,ali,treefile):
        
        self.treefile  = treefile
        self.alignment = ali
        if self.alignment == None:
            raise AttributeError('No alignment defined.')
        if not os.path.isfile(self.treefile):
            raise IOError('Tree file was not found.')
    
    def __inst__(self):
        
        return 'echo "y" | %s %s %s -wsl' % (
            TREEPUZZLE_BIN,self.alignment.getPhylip(),self.treefile)
    
    def getSiteLikelihoodFile(self):
        
        sproc = system(self.__inst__(
            ),stdout=PIPE,stderr=PIPE,shell=True)
        self._output = sproc.communicate()[0]
        if not os.path.isfile('%s.sitelh' % (self.treefile)):
            raise IOError('TREE-PUZZLE did not create site-likelihood output.')
        else: return '%s.sitelh' % (self.treefile)
    
class consel(object): # TODO: Move over to executable interface.
    
    ''' Wrap CONSEL for confidence analysis of a set of trees comprising a
    phylogenetic landscape (landscape.landscape object). '''
    
    def __init__(self,landscape,name):
        
        self.landscape = landscape
        self.name      = name
        self.alignment = landscape.alignment
        self.sitelh    = None
        self.raw       = None
        self.rmt       = None
        self.pv        = None
        self.auvals    = list()
        self.interval  = list()
        self.rejected  = list()
    
    def _convertRawData(self):
        
        if not os.path.isfile('%s.mt' % (self.name)):
            sproc = system('seqmt --puzzle %s %s.mt' % (self.sitelh,self.name),
                           stdout=PIPE,stderr=PIPE,shell=True)
            self._out = sproc.communicate()[0]
            if not os.path.isfile('%s.mt' % (self.name)):
                raise IOError('CONSEL seqmt did not create mt file.')
        self.raw = '%s.mt' % (self.name)
    
    def _createReplicates(self):
        
        if not os.path.isfile('%s.rmt' % (self.name)):
            sproc = system('makermt %s' % (self.raw),stdout=PIPE,
                           stderr=PIPE,shell=True)
            self._out = sproc.communicate()[0]
            if not os.path.isfile('%s.rmt' % (self.name)):
                raise IOError('CONSEL markermt did not create rmt file.')
        self.rmt = '%s.rmt' % (self.name)
    
    def _run(self):
        
        if not os.path.isfile('%s.pv' % (self.name)):
            sproc = system('consel %s.rmt' % (self.name),stdout=PIPE,
                           stderr=PIPE,shell=True)
            self._out = sproc.communicate()[0]
        self.pv = '%s.pv' % (self.name)
        
    def _getAU(self):
        
        if not os.path.isfile(self.pv):
            raise IOError('CONSEL consel p-value file not found.')
        
        sproc = system('catpv %s' % (self.pv),stdout=PIPE,
                       stderr=PIPE,shell=True)
        pvout = sproc.communicate()[0]
        for line in pvout.split('\n'):
            spl = line.split()
            if len(spl) < 3: continue
            elif not spl[2].isdigit(): continue
            it = (int(spl[2])-1,float(spl[4]))
            self.auvals.append(it)
            if it[1] >= 0.05: self.interval.append(
                self.landscape.getNodeNames()[it[0]])
            else: self.rejected.append(
                self.landscape.getNodeNames()[it[0]])
        self._out = pvout
            
    def compute(self):
        
        ''' Compute the AU test confidence interval for a landscape. '''
        
        with aTemporaryDirectory():
            if not os.path.isfile('%s.trees.sitelh' % (self.name)):
                self.treefile  = self.landscape.toTreeFile(self.name + '.trees',False)
                self.treepuzz  = treepuzzle(self.alignment,self.treefile)
                self.sitelh = self.treepuzz.getSiteLikelihoodFile()
            else: self.sitelh = '%s.trees.sitelh' % (self.name)
            self._convertRawData()
            self._createReplicates()
            self._run()
            self._getAU()