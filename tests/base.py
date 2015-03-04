from pylogeny.newick import shuffleLeaves, newickParser
from pylogeny.tree import tree as treeObject
from pylogeny.alignment import phylipFriendlyAlignment as alignment
from pylogeny.landscape import landscape
from pylogeny.landscapeWriter import landscapeWriter, landscapeParser
from unittest import TestLoader as loader, TestCase as testCase, TextTestRunner as tests
from os.path import isfile
from os import unlink

TESTS_ALIGNMENT = 'al.fasta'
TESTS_OPERATOR  = 'SPR'

class phylogeneticLandscapeTest(testCase):
    
    landscape = None
    
    @classmethod
    def setUpClass(cls):
        if cls.landscape == None:
            ali = alignment(TESTS_ALIGNMENT)
            cls.landscape = landscape(ali,operator=TESTS_OPERATOR)
            
    @classmethod
    def tearDownClass(cls):
        if isfile('al.landscape'): unlink('al.landscape')    
    
    def assertHasNoDuplicateTrees(self):                                
        for i in self.landscape.iterNodes():             
            t = self.landscape.getTree(i)
            f = self.landscape.findTreeTopologyByStructure(t.getStructure())
            if (f == None): continue
            self.assertFalse((f != i and t == self.landscape.getTree(f)),
                             msg="Tree objects (%s,%s) are equal but have different IDs." % (str(f),str()))
            self.assertFalse((f != i),msg="Tree %s duplicate to %s." % (str(f),str(i)))
    
    def resetLandscape(self):    
        ali = alignment(TESTS_ALIGNMENT)
        self.landscape = landscape(ali,operator=TESTS_OPERATOR)
    
    def getRandomTree(self,treeid):
        treeNewick = self.landscape.getTree(treeid).getNewick()
        parserObj  = newickParser(treeNewick)
        shuffleLeaves(parserObj.parse())
        return treeObject(str(parserObj.parsed_structure)+';',check=True)
    
    def runHeuristic(self,heu,start=0,*args):
        h = heu(self.landscape,self.landscape.getNode(start),*args)
        h.explore()
        self.assertGreater(len(self.landscape),1)
        self.assertHasNoDuplicateTrees()    
