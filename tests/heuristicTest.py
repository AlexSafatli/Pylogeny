from base import *
from Odin import heuristic

class heuristicTest(phylogeneticLandscapeTest):

    def runLinearHeuristic(self,heu,start=0):
        if issubclass(heu,heuristic.phylogeneticLinearHeuristic):
            self.runHeuristic(heu,start)
            return
        self.fail('Not a proper heuristic.')

class implementedHeuristicsTest(heuristicTest):

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

    suite = loader().loadTestsFromTestCase(implementedHeuristicsTest)
    tests(verbosity=2).run(suite)
