from unittest import TestLoader as loader, TestCase as testCase, TextTestRunner as tests
from os.path import isfile
from pylogeny.tree import tree, treeSet
from pylogeny.alignment import alignment

class treeTest(testCase):

    alignm = None
    tree_  = None
    trees  = list()
    
    @classmethod
    def setUpClass(cls):
        if cls.alignm == None:
            cls.alignm = alignment('al.fasta')
            cls.tree_  = cls.alignm.getApproxMLTree()
            cls.trees  = treeSet()
            cls.trees.addTree(cls.tree_)
    
    @classmethod
    def tearDownClass(cls): pass
    
    def test_tree_init(self):
        
        t = tree(self.tree_.getNewick(),check=True)
        self.assertTrue(type(t) == tree)
    
    def test_tree_initNoCheck(self):
        
        t = tree(self.tree_.getNewick(),check=False)
        self.assertTrue(type(t) == tree)
        self.assertEqual(t.newick,self.tree_.getNewick())
    
    def test_tree_toTopology(self):
         
        t = self.tree_.toTopology()
        self.assertEqual(t.toNewick(),self.tree_.getNewick())
        
    def test_treeSet_indexOfTree(self):
        
        self.assertGreaterEqual(self.trees.indexOf(self.tree_),0)
        
    def test_treeSet_addTree(self):
        
        l = len(self.trees)
        t = self.tree_
        self.trees.addTree(t)
        self.assertEqual(len(self.trees),l+1)
    
    def test_treeSet_removeTree(self):
        
        l = len(self.trees)
        t = self.tree_
        self.trees.removeTree(t)
        self.assertEqual(len(self.trees),l-1)
        
if __name__ == '__main__':
    suite = loader().loadTestsFromTestCase(treeTest)
    tests(verbosity=2).run(suite)
