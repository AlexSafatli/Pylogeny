''' Defines an interface to manage interfacing with the system for respective application calls and implements multiple of these for executables such as FastTree and RAxML. Requires a UNIX environment. '''

# Date:   Oct 16 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import tree
from abc import ABCMeta as abstractclass, abstractmethod
from os import mkdir, getcwd, chdir
from os.path import abspath, isdir, isfile
from subprocess import call, PIPE, Popen as system

# Constants

E_FASTTREE = 'fasttree'
E_RAXML    = 'raxmlHPC'
E_TREEPUZZ = 'puzzle'

# Executable Existence Function

def exeExists(cmd):
    
    ''' Determines whether a function exists in a UNIX environment. '''
    
    return call('type %s'%(cmd),shell=True,stdout=PIPE,stderr=PIPE) == 0

# Temporary Directory Context

class aTemporaryDirectory(object):
    
    ''' A class intended to be used as a context manager that allows
    Python to run in a temporary directory for a finite period of time. '''
    
    def __init__(self,dir=None):
        
        self.path = '/tmp' + dir
        if not isdir(self.path): mkdir(self.path)
        self.current = getcwd()
        
    def __enter__(self):      chdir(self.path)
    def __exit__(self,*args): chdir(self.current)

# Executable Interface

class executable(object):
    
    ''' An interface for the instantation and running of a single instance for a given application. '''

    exeName = None
    __metaclass__ = abstractclass
    
    @abstractmethod
    def getInstructionString(self): pass
    
    def run(self):
        
        ''' Perform a run of this application. '''
        
        if (self.exeName and not exeExists(self.exeName)):
            raise SystemError('%s is not installed on your system.'%(self.exeName))
        sproc = system(self.getInstructionString(),
                       stdout=PIPE,stderr=PIPE,shell=True)
        return sproc.communicate()[0]

# Executable Classes

class treepuzzle(executable):
    
    ''' Wrap TREE-PUZZLE in order to create an intermediate file for CONSEL to read and assign confidence to a set of trees. Requires TREE-PUZZLE to be installed. '''

    exeName = E_TREEPUZZ

    def __init__(self,ali,treefile):
        self.treefile  = treefile
        self.alignment = ali
        if self.alignment == None:
            raise AttributeError('No alignment defined.')
        
    def getInstructionString(self):
        return 'echo "y" | %s %s %s -wsl' % (
            self.exeName,self.alignment.getPhylip(),self.treefile)
    
    def getSiteLikelihoodFile(self):
        self._output = self.run()
        if not isfile('%s.sitelh' % (self.treefile)):
            raise IOError('TREE-PUZZLE did not create site-likelihood output.')
        else: return '%s.sitelh' % (self.treefile) 

class consel(executable):
    
    ''' Denotes a single run of the CONSEL workflow in order to acquire a confidence interval and perform an AU test on a set of trees. Requires CONSEL to be installed. '''
    
    def __init__(self,treeset,alignment,name):
        
        self.treeset   = treeset
        self.name      = name
        self.alignment = alignment
        self.sitelh    = None
        self.raw       = None
        self.rmt       = None
        self.pv        = None
        self.auvals    = list()
        self.interval  = tree.treeSet()
        self.rejected  = tree.treeSet()
        self.instruction = ''
        
    def getInstructionString(self): return self.instruction
    
    def _convertRawData(self):
        
        if not isfile('%s.mt' % (self.name)):
            self.instruction = 'seqmt --puzzle %s %s.mt' % (self.sitelh,self.name)
            self._out = self.run()
            if not isfile('%s.mt' % (self.name)):
                raise IOError('CONSEL seqmt did not create mt file.') 
        self.raw = '%s.mt' % (self.name)
        
    def _createReplicates(self):
        
        if not isfile('%s.rmt' % (self.name)):
            self.instruction = 'makermt %s' % (self.raw)
            self._out = self.run()
            if not isfile('%s.rmt' % (self.name)):
                raise IOError('CONSEL markermt did not create rmt file.')
        self.rmt = '%s.rmt' % (self.name)
        
    def _run(self):
        
        if not isfile('%s.pv' % (self.name)):
            self.instruction = 'consel %s.rmt' % (self.name)
            self._out = self.run()
        self.pv = '%s.pv' % (self.name)
        
    def _getAU(self):
        
        self.instruction = 'catpv %s' % (self.pv)
        pvout = self.out()
        for line in pvout.split('\n'):
            spl = line.split()
            if len(spl) < 3: continue
            elif not spl[2].isdigit(): continue
            it = (int(spl[2])-1,float(spl[4]))
            self.auvals.append(it)
            if it[1] >= 0.05: self.interval.addTree(self.treeset[it[0]])
            else: self.rejected.addTree(self.treeset[it[0]])    
        self._out = pvout
    
    def getInterval(self):
        
        ''' Compute the AU test. Return the interval of trees. '''
        
        with aTemporaryDirectory('/consel__py/'):
            if not isfile('%s.trees.sitelh' % (self.name)):
                self.treefile = self.treeset.toTreeFile(self.name + '.trees',False)
                self.treepuzz = treepuzzle(self.alignment,self.treefile)
                self.sitelh   = self.treepuzz.getSiteLikelihoodFile()
            else: self.sitelh = '%s.trees.sitelh' % (self.name)
            self._convertRawData()
            self._createReplicates()
            self._run()
            self._getAU()
        return self.interval

class fasttree(executable):
    
    ''' Denotes a single run of the FastTree executable in order to acquire an approximate maximum likelihood tree for the input alignment. See http://www.microbesonline.org/fasttree/ for more information on FastTree. Requires FastTree to be installed. '''
    
    exeName = E_FASTTREE
    
    def __init__(self,inp_align,out_file=None,isProtein=True):
        self.alignment = inp_align
        self.isProtein = isProtein
        self.out       = out_file
        
    def getInstructionString(self):
        s = '%s -quiet -nopr < %s' % (self.exeName,self.alignment)
        if self.out: s += ' > %s' % (self.out)
        return s
    
class raxml(executable):
    
    ''' Denotes a single run of the RAxML executable. See http://sco.h-its.org/exelixis/software.html for more information on RAxML. Requires RAxML to be installed. '''
    
    exeName = E_RAXML
    
    def __init__(self,inp_align,out_file,model=None,
                 is_Protein=True,interTrees=False,
                 alg=None,startingTree=None,rapid=False,
                 slow=False,optimizeBootstrap=False,numboot=100,
                 log=None,wdir=None):
        self.alignment     = inp_align
        self.out           = out_file
        self.model         = model
        self.protein       = is_Protein
        self.intermediates = interTrees
        self.alg           = alg
        self.startingTree  = startingTree
        self.slowBoot      = slow
        self.rapidBoot     = rapid
        self.numBoot       = numboot
        self.optBoot       = optimizeBootstrap
        self.log           = log
        self.workdir       = abspath(wdir)
        if not self.model:
            if is_Protein: self.model = 'PROTGAMMAJTT'
            else:          self.model = 'GTRGAMMA'
            
    def getInstructionString(self):
        s = '%s -s %s -n %s -m %s' % (self.exeName,self.alignment,self.out,self.model)
        if self.alg: s += ' -f %s' % (self.alg)
        if self.intermediates: s += ' -j'
        if self.startingTree: s += ' -t %s' % (self.startingTree)
        if self.rapidBoot or self.slowBoot:
            if self.rapidBoot:  s += ' -x 234534251 -N %d' % (self.numBoot)
            elif self.slowBoot: s += ' -b 234534251 -N %d' % (self.numBoot)
        if self.optBoot: s += ' -k'
        if self.workdir: s += ' -w %s' % (self.workdir)
        if self.log: s += ' > %s' % (self.log)
        return s
    
    def runFunction(self,alg):
        self.alg = alg
        self.run()

