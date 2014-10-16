''' Python interface for the FastTree executable. See http://www.microbesonline.org/fasttree/ for more information on FastTree. Requires FastTree to be installed. '''

# Date:   Jan 27 2013
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from subprocess import PIPE, Popen as system
FT_EXECUTABLE = 'fasttree'

def exeExists(cmd):
    return subprocess.call('type %s' % (cmd),shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) == 0

class fasttree:
    
    ''' Denotes a single run of the FastTree executable. '''
    
    def __init__(self,inp_align,out_file=None,isProtein=True):
        self.alignment = inp_align
        self.isProtein = isProtein
        self.out       = out_file
        
    def getInstructionString(self):
        s = '%s -quiet -nopr < %s' % (FT_EXECUTABLE,self.alignment)
        if self.out: s += ' > %s' % (self.out)
        return s
    
    def run(self):
        
        ''' Perform a run of the FastTree executable in order to acquire
        an approximate maximum likelihood tree for the input alignment. '''
        
        if not exeExists(FT_EXECUTABLE):
            raise SystemError('FastTree is not installed on your system.')
        sproc = system(self.getInstructionString(),
                       stdout=PIPE,stderr=PIPE,shell=True)
        if not self.out: return sproc.communicate()[0]
        else:            return self.out
        