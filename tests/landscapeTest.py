from base import *
from random import sample

class landscapeTest(phylogeneticLandscapeTest):

    def test_findTreeTopologyByStructure(self):
        tree = self.landscape.getTree(0)
        find = self.landscape.findTreeTopologyByStructure(tree.getStructure())
        self.assertFalse(find == None)
        self.assertTrue(type(self.landscape.getTree(find)) == treeObject)
        self.assertEquals(0,find)

    def _exploreTreeByName(self,treeid):
        self.landscape.exploreTree(treeid)
        self.assertGreater(len(self.landscape),1)
        for tr in self.landscape.getNeighborsFor(treeid):
            self.assertIn(tr,self.landscape.getNodeNames())
            t = self.landscape.getTree(tr)
            self.assertIsInstance(t,treeObject)
            self.assertEqual(self.landscape.findTreeTopologyByStructure(t.getStructure()),tr)
            self.assertIsNotNone(self.landscape.getEdge(treeid,tr))        
        
    def test_getTreeNewick(self):
        treeNewick = self.landscape.getTree(0).getNewick()
        self.assertTrue(type(treeNewick) == str and treeNewick != None)
        self.assertTrue(len(treeNewick) > 0)
    
    def test_addDuplicateTree(self):
        randomTree = None
        while (randomTree == None or 
               self.landscape.findTreeTopologyByStructure(randomTree.getStructure()) != None):
            randomTree = self.getRandomTree(0)
            dupTree    = treeObject(randomTree.getStructure())
        self.landscape.addTree(randomTree)
        self.assertIsNotNone(self.landscape.findTreeTopologyByStructure(randomTree.getStructure()))
        self.assertRaises(AssertionError,self.landscape.addTree,dupTree)
    
    def test_addRemoveTree(self):
        randomTree = self.getRandomTree(0)
        self.landscape.addTree(randomTree)
        self.assertIsNotNone(self.landscape.findTreeTopologyByStructure(randomTree.getStructure()))
        self.landscape.removeTree(randomTree)
        self.assertIsNone(self.landscape.findTreeTopologyByStructure(randomTree.getStructure()))

    def test_exploreRootTree(self):
        self._exploreTreeByName(0)
    
    def test_exploreRandomTree(self):
        randomTree = self.getRandomTree(0)
        i = self.landscape.addTree(randomTree)
        self._exploreTreeByName(i)
        
    def test_getBipartitions(self):
        for node in sample(self.landscape.getNodeNames(),min([len(self.landscape),10])):
            bps = self.landscape.getVertex(node).getBipartitions()
            self.assertGreater(len(bps),0)
    
    def test_getProperTrees(self):
        stringList = []
        for i in self.landscape.iterNodes():
            stringList.append(self.landscape.getVertex(i).getProperNewick())
        self.assertEquals(len(stringList),len(self.landscape))

    def test_getPathOfBestImprovement(self):
        ids = [0] + self.landscape.getPathOfBestImprovement(0)
        ml = lambda i: self.landscape.getTree(i).getScore()[0]
        self.assertTrue(all([ml(ids[i]) <= ml(ids[i+1]) for i in xrange(1,len(ids))]))
    
    def test_readAndWrite(self):
        struct = lambda i,l: l.getTree(i).getStructure()
        writer = landscapeWriter(self.landscape,'al')
        writer.writeFile()
        self.assertTrue(isfile('al.landscape'))
        reader = landscapeParser('al.landscape')
        l,n = reader.parse()
        self.assertEqual(len(l),len(self.landscape))
        for i in self.landscape.getNodeNames():
            self.assertEquals(struct(i,l),struct(i,self.landscape))

if __name__ == '__main__':

    suite = loader().loadTestsFromTestCase(landscapeTest)
    tests(verbosity=2).run(suite)
