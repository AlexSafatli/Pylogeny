# Patricia Tree Testing

import platform, random
from unittest import TestLoader as loader, TestCase as testCase, TextTestRunner as tests
from base import phylogeneticLandscapeTest
from pylogeny.landscape import landscape
from pylogeny.base import patriciaTree
from pylogeny import tree as otree
from pylogeny.landscapeWriter import landscapeWriter, landscapeParser
from pylogeny.alignment import phylipFriendlyAlignment as alignment
from pylogeny.heuristic import parsimonyGreedy

NUM_WORDS_TO_TEST = 60000

class patriciaTreeTest(phylogeneticLandscapeTest):

    words       = None
    wordPatTree = None
    newiPatTree = None
    
    @classmethod
    def setUpClass(cls):
        
        # Acquire words from the unix dictionary.
        o = open('/usr/share/dict/words')
        raw_words = o.readlines()
        o.close()
        cls.words = list()
        for num in xrange(NUM_WORDS_TO_TEST): # Try a finite number of distinct words.
            word = random.choice(raw_words).strip('\n')
            while (word in cls.words):
                word = random.choice(raw_words).strip('\n')
            cls.words.append(word)
        
        # Set up the landscape for another source of data.
        ali = alignment('tests/al.fasta')
        cls.landscape = landscape(ali,starting_tree=ali.getApproxMLTree(),
                                   root=True,operator='SPR')
        
        # Set up two PATRICIA trees.
        cls.wordPatTree = patriciaTree()
        cls.newiPatTree = patriciaTree()
    
    @classmethod
    def tearDownClass(cls):
        
        # Write tree string representations to files.
        wordPatString = str(cls.wordPatTree)
        newiPatString = str(cls.newiPatTree)
        o = open('wordPatTree','w')
        o.write(wordPatString)
        o.close()
        o = open('newiPatTree','w')
        o.write(newiPatString)
        o.close()
       
        # Set all variables to None.
        cls.words = None
        cls.landscape = None
        cls.wordPatTree = None
        cls.newiPatTree = None
    
    def test_addSingleWordString(self):
        self.wordPatTree.insert('Door')
        self.assertTrue('Door' in self.wordPatTree)
    
    def test_addMultipleWordStrings(self):
        prevlen = len(self.wordPatTree)
        for word in self.words:
            currlen = len(self.wordPatTree)
            self.wordPatTree.insert(word)
            self.assertTrue(word in self.wordPatTree)
            self.assertEquals(len(self.wordPatTree),currlen+1)
        self.assertEquals(len(self.wordPatTree),prevlen+len(self.words))
        
    def test_addSingleNewickTree(self):
        newick = self.landscape.getTree(0).toNewick()
        self.newiPatTree.insert(newick)
        self.assertTrue(newick in self.newiPatTree)
    
    def test_addMultipleNewickTrees(self):
        self.runHeuristic(parsimonyGreedy)
        prevlen = len(self.newiPatTree)
        for tree in self.landscape.iterTrees():
            self.assertTrue(type(tree) == otree.tree)
            if tree == self.landscape.getRootTree(): continue
            currlen = len(self.newiPatTree)
            self.newiPatTree.insert(tree.toNewick())
            self.assertTrue(tree.toNewick() in self.newiPatTree)
            self.assertEquals(len(self.newiPatTree),currlen+1)
        self.assertEquals(len(self.newiPatTree),prevlen+len(self.landscape)-1)

if __name__ == '__main__':
    # Ensure on a UNIX platform. Naively using OS word list from file.
    os = platform.system()
    if os != 'Linux' and os != 'Darwin': # (Linux or Mac OS X)
        raise OSError('Not on an applicable platform! Current platform: %s' % (os))
        
    # Set up tests.
    suite = loader().loadTestsFromTestCase(patriciaTreeTest)
    tests(verbosity=2).run(suite)
