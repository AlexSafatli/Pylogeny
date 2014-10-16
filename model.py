''' Phylogenetic tree scoring models; intended to be coupled with the use of pytbeaglehon (BEAGLE) high-performance library. '''

# Date:   Jun 9 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

# High-Performance Numeric Computation Library (pandas)

import pandas as P

# pytbeaglehon coupling

pytbeaglehonEnabled = True
try:
    import pytbeaglehon.disc_state_cont_time_model as BEAGLEModels
    from pytbeaglehon.disc_char_type import DNAType, AAType
    from pytbeaglehon.tests.util import TreeForTesting as BEAGLETree
except: pytbeaglehonEnabled = False

# Exception Handling

class PhyloModelError(Exception):
    def __init__(self,v): self.value = v
    def __str__(self): return repr(self.value)

# Model a biological sequence as a discrete state model.

class DiscreteStateModel(object):
    
    ''' Initialize a discrete state model 
    for phylogenetic data. State frequencies and character
    time are determined from the given alignment object. '''
    
    def __init__(self,alignment):
        
        self.alignment    = alignment
        self.charType     = None
        self.stateList    = None
        self.seqMatrix    = None
        self.rawFreqMap   = None
        self.rawFreqs     = None
        self.stateFreqs   = None
        self._totalFreq   = 0.0
        if pytbeaglehonEnabled:
            self._determineCharType()
            self._transformAlignmentToMatrix()
            self._determineStateFreqs()
        
    # Internals
    
    def _determineCharType(self):
        
        ''' PRIVATE: Determine the character type intrinsic to the
        biological sequence alignment given to the instance. '''
        
        a = self.alignment
        if (a.getDim() == 4): self.charType = DNAType()
        elif (a.getDim() >= 20 and a.getDim() <= 21):
            self.charType = AAType()
        else: raise PhyloModelError(
            'Unrecognizable character type in alignment.')
    
    def _transformAlignmentToMatrix(self):
        
        ''' PRIVATE: Transform instance alignment to a matrix of 
        states as integers. '''
        
        lioflists = []
        strList   = map(str.upper,self.alignment.toStrList())
        translate = self.charType.to_indices
        for l in strList: lioflists.append(translate(l))
        self.stateList = lioflists
        self.seqMatrix = P.DataFrame(lioflists)
    
    def _determineStateFreqs(self):
        
        ''' PRIVATE: Compute the state frequencies for all states if
        character type is DNA. Otherwise, leave undefined. '''
    
        if self.charType.get_num_states() != 4: return # DNAType?
        self.rawFreqMap = {}
        
        # Get all frequencies from the sequence matrix dataframe.
        freqS = self.seqMatrix.stack().value_counts()
        
        # Get all present states that are found in the dataframe.
        present = freqS.keys()
        
        # Get all possible states.
        states = self.charType.get_states()
        
        # Enumerate over all possible states.
        for state in states:
            i = self.charType.to_indices([state])[0]
            if i in present:
                self.rawFreqMap[state] = freqS.get(i)
            else: self.rawFreqMap[state] = 0
        self.rawFreqs = self.rawFreqMap.items() # Flatten.
        
        # Calculate actual state frequencies.
        total = float(sum([i[1] for i in self.rawFreqs]))
        self._totalFreq = total
        self.stateFreqs = map(lambda x:x[1]/total,self.rawFreqs)
    
    # Accessors    
    
    def getAlignment(self): return self.alignment
    def getAlignmentAsStateList(self): return self.stateList
    def getSequenceMatrix(self): return self.seqMatrix
    def getCharType(self): return self.charType
    def getStateFreqs(self): return self.stateFreqs
    def getRawStateFreqs(self): return self.rawFreqs
    def getRawStateFreqsAsList(self): return self.getRawStateFreqs()
    def getRawStateFreqsAsDict(self): return self.rawFreqMap
    
    def getFrequencyOfState(self,i):
        if self.rawFreqMap:
            if i in self.rawFreqMap:
                return self.rawFreqMap[i]/self._totalFreq
    
    def getRawFrequencyOfState(self,i):
        if self.rawFreqMap:
            if i in self.rawFreqMap:
                return self.rawFreqMap[i]