''' Serialize a phylogenetic landscape into a JSON object. '''

# Date:   Jan 27 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from landscapeWriter import landscapeWriter
import json

class JSONWriter(landscapeWriter):
    
    def __init__(self, ls, name):
        super(JSONWriter,self).__init__(ls,name)

    def nodeToJSON(self,node):

        n = self.landscape.getVertex(node)
        ml,pars = n.getScore()
        og      = n.getOrigin()
        lmax    = n.isLocalOptimum()
        mimpr   = n.getBestImprovement()
        expl    = n.isExplored()
        dg      = len(n.getNeighbors())
        val     = int(n.isFailed())
        return {'name':str(node),'group':int(node),
               'ml':ml,'parsimony':pars,'origin':og,
               'explored':expl,'degree':dg,
               'isLocalMax':lmax,'maximprov':str(mimpr),
               'value':val} 

    def _jsonNodes(self):
        
        for node in self.graph.node:
            if self.graph.degree(node) > 0:
                yield self.nodeToJSON(node)
    
    def _jsonEdges(self,jsonNodes):
        
        nodeli = [x['group'] for x in jsonNodes]
        for node in nodeli:
            neighbors = self.graph.neighbors(node)
            for neigh in neighbors:
                if neigh in nodeli:
                    yield {'source':nodeli.index(node),
                           'target':nodeli.index(neigh),
                           'value': 1}    
    
    def _jsonGraph(self,nodes=None):
        
        # Initialize variables.
        if not nodes: nodes = self._jsonNodes()
        ls       = self.landscape
        locks    = ls.locks
        toVis    = []
        allScore = True

        # Check scores.
        for node in nodes:

            # Get node's log-likelihood score.
            nS = node['ml']
            if (nS == None): allScore = False
            else: toVis.append(node)

        if len(toVis) == 0: allScore = False        
        return toVis, allScore
    
    def getOnlyImprovements(self,groups=None):
        
        # Are there any groups to uncollapse?
        if (groups == None): groups = list()
        
        # Get all relevant nodes.
        nodes,inds = [self.nodeToJSON(0)],[0]
        paths = self.landscape.iterAllPathsOfBestImprovement()
        for path in paths:
            for it in path:
                if not it in inds:
                    nodes.append(self.nodeToJSON(it))
                    inds.append(it)
        
        # Get node and link information.
        nodes, aS = self._jsonGraph(nodes)
        links = [x for x in self._jsonEdges(nodes)]
        out   = json.dumps(
            {"nodes":nodes,"links":links,"graphid":self.name,
            "scored":aS},sort_keys=True,indent=4,separators=(',',': '))
        return out
    
    
    def getCompleteLandscape(self): return self.getJSON()

    def getJSON(self):

        # Get node and link information.
        nodes, aS = self._jsonGraph()
        nodes = [node for node in nodes if node['ml'] != None]
        links = [x for x in self._jsonEdges(nodes)]
        
        # Output to JSON string.
        out = json.dumps({"nodes":nodes,"links":links,"graphid":self.name,
                         "scored":aS},sort_keys=True,indent=4,separators=(',',': '))
        return out
        