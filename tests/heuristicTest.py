from base import *
from Odin import heuristic

class heuristicTest(pylogenyTest):

    landscape = None
    
    def runLinearHeuristic(self,heu,start=0):
        if issubclass(heu,heuristic.phylogeneticLinearHeuristic):
            h = heu(self.landscape,self.landscape.getNode(start))
            h.explore()
            self.assertGreater(len(self.landscape),1)
            return
        self.fail('Not a proper heuristic.')

    def test_runParsimonyGreedyHeuristic(self):
        self.resetLandscape()
        self.runLinearHeuristic(heuristic.parsimonyGreedy)
    
    def test_runLikelihoodGreedyHeuristic(self):
        self.resetLandscape()
        self.runLinearHeuristic(heuristic.likelihoodGreedy)
        
    def test_runSmoothGreedyHeuristic(self):
        self.resetLandscape()
        self.runLinearHeuristic(heuristic.smoothGreedy)

if __name__ == '__main__':

    suite = loader().loadTestsFromTestCase(heuristicTest)
    tests(verbosity=2).run(suite)
