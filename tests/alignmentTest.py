from base import *
from p4 import Sequence

class alignmentTest(testCase):
    
    alignment = None
    
    @classmethod
    def setUpClass(cls):
        if cls.alignment == None:
            cls.alignment = alignment(TESTS_ALIGNMENT)
    
    def test_toString(self):
        self.assertIsInstance(str(self.alignment),str)
        lines = str(self.alignment).split('\n')
        for i in xrange(len(lines)):
            self.assertEqual(lines[i],self.alignment.getSequence(i))    
    
    def test_getLength(self):
        self.assertEquals(len(self.alignment),len(self.alignment.data))
        self.assertEquals(len(self.alignment),self.alignment.getSize())
    
    def test_getItem(self):
        for seq in self.alignment:
            self.assertIsNotNone(seq)
            self.assertTrue(type(seq) == Sequence)
            
if __name__ == '__main__':

    suite = loader().loadTestsFromTestCase(alignmentTest)
    tests(verbosity=2).run(suite)