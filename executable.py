''' Defines an interface to manage interfacing with the system for respective application calls and implements multiple of these for executables such as FastTree and RAxML. '''

# Date:   Oct 16 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from os.path import abspath
from subprocess import PIPE, Popen as system

# Constants

E_FASTTREE = 'fasttree'
E_RAXML    = 'raxmlHPC'

# Executable Existence Function

def exeExists(cmd):
    return subprocess.call('type %s'%(cmd),shell=True,stdout=PIPE,stderr=PIPE) == 0

# Executable Interface

class executable(object):
    
    ''' An interface for the instantation and running of a single instance for a given application. '''

    exeName = None
    
    def __init__(self): pass
    
    def getInstructionString(self): return ''
    
    def run(self):
        
        ''' Perform a run of this application. '''
        
        if not exeExists(self.exeName):
            raise SystemError('%s is not installed on your system.'%(self.exeName))
        sproc = system(self.getInstructionString(),
                       stdout=PIPE,stderr=PIPE,shell=True)
        if not self.out: return sproc.communicate()[0]
        else:            return self.out    

# Executable Classes

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
        