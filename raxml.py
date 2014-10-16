''' Python interface for RAxML executable. '''

# Date:   Jan 27 2013 
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from os.path import abspath
from subprocess import PIPE, Popen as system

RX_EXECUTABLE = 'raxmlHPC'

class raxml:
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
        s = '%s -s %s -n %s -m %s' % (RX_EXECUTABLE,self.alignment,self.out,self.model)
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
    def run(self):
        sproc = system(self.getInstructionString(),
                       stdout=PIPE,stderr=PIPE,shell=True)
        return sproc.communicate()[0]
        