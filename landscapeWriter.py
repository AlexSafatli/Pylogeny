''' Serialize a phylogenetic landscape into a Python-readable file (Pickle). '''

# Date:   Apr 9 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import os, pickle
from landscape import landscape

class landscapeWriter(object):

    ''' Encapsulate the writing of a landscape
    to a file format. '''

    def __init__(self, landscape, name):
        
        # Fields
        self.landscape = landscape
        self.name      = name
        self.cleankey  = ''
        self.cleanpath = ''
        self.graph     = landscape.graph
        
        # Clean the input landscape name.
        self._cleanName()

    def _cleanName(self):

        for x in self.name:
            if x.isalnum(): self.cleankey += x
            else: self.cleankey += '_'
        self.cleanpath = '%s.landscape' % (self.cleankey)

    def _fixDataStructures(self):
        
        # Change the list proxy found in the data structure to a list.
        self.landscape._locks = self.landscape.locks
        self.landscape.locks = [x for x in self.landscape.locks]

        # Remove any temporary files from the alignment object.
        if self.landscape.alignment:
            self.landscape.alignment.close()
            if hasattr(self.landscape.alignment,'_pllmodel'):
                del self.landscape.alignment._pllmodel

    def _resetDataStructures(self):
        
        # Return structures to original forms.
        self.landscape.locks = self.landscape._locks
        if self.landscape.alignment:
            self.landscape.alignment = self.landscape.alignment.recreateObject()

    def _dump(self,path='.'):

        # Fix structures in landscape.
        self._fixDataStructures()

        # Open the file.
        fpath = os.path.join(path,self.cleanpath)
        o = open(fpath, mode='wb')
        
        # Dump the landscape.
        pickle.dump((self.landscape,self.name),o)
        
        # Close the file.
        o.close()
        
        # Return structures to previous form.
        self._resetDataStructures()
        
        # Return the path.
        return fpath

    def writeFile(self,path='.'):
        
        ''' Write the landscape serialized file to given path. '''

        return self._dump(path)

class landscapeParser:

    ''' Encapsulates the construction of a landscape
    object from a pickle file. '''

    def __init__(self,path):
        self.file = path

    def _onparse(self,landsc):
        if landsc.alignment: landsc.alignment.paths = {}         

    def parse(self):
        
        ''' Parse the file. '''

        # Get landscape object ready.
        landsc = None

        # Read the file.
        o = open(self.file,'rb')
        landsc,name = pickle.load(o)
        self._onparse(landsc) # Do any fixes.
        o.close()

        return (landsc,name)