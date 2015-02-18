''' Serialize a phylogenetic landscape into an SQLlite database file made up of
three components: all tree IDs and respective scores, the alignment file as a
set of sequences, and a representation of the graph as an edge list. '''

# Date:   Apr 9 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from __version__ import VERSION
from alignment import phylipFriendlyAlignment as alignment
from landscape import landscape
from database import SQLiteDatabase
import os, tree

class landscapeWriter(object):

    ''' Encapsulate the writing of a landscape to a file format. '''

    def __init__(self, landscape, name):
        
        # Fields
        self.landscape = landscape
        self.name      = name
        self.cleankey  = ''
        self.graph     = landscape.graph
        self.database  = None        
        
        # Clean the input landscape name.
        self._cleanName()

    def _cleanName(self):

        for x in self.name:
            if x.isalnum(): self.cleankey += x
            else: self.cleankey += '_'
        self.cleansuff = '%s.landscape' % (self.cleankey)

    def _schema(self):
        
        # Get the database object.
        dbobj = self.database
        
        # Create all content tables.
        dbobj.newTable('alignment',('seqid','integer'),('seqname','text'),('sequence','text'))
        dbobj.newTable('trees',('treeid','integer'),('name','text'),('newick','text'),
                       ('origin','text'),('ml','real'),('pars','real'),('explored','boolean'))
        dbobj.newTable('graph',('source','integer'),('origin','integer'))
        dbobj.newTable('locks',('treeid','integer'),('branchid','integer'))
        
        # Create a table to store metadata.
        dbobj.newTable('metadata',('key','text'),('value','text'))

    def _dump(self,path='.'):

        # Open the file.
        fpath = os.path.join(path,self.cleansuff)
        if os.path.isfile(fpath): os.unlink(fpath)
        o = SQLiteDatabase(fpath)
        self.database = o
        
        # Construct the schema for the landscape.
        self._schema()
        
        # Include metadata about the landscape.
        o.insertRecord('metadata',['name',self.name])
        o.insertRecord('metadata',['version',VERSION])
        
        # Add the alignment.
        index = 0
        if self.landscape.alignment != None:
            for s in self.landscape.alignment:
                o.insertRecord('alignment',[index,s.name,s.sequence])
                index += 1
        
        # Add all of the trees.
        for i in self.landscape.iterNodes():
            t = self.landscape.getTree(i)
            find = self.landscape.findTreeTopologyByStructure(t.getStructure())
            if (find == None):
                raise AssertionError('Tree %s was in landscape but no stored structure.' % (str(i)))
            elif (find == i):
                # Verify unique identity of this tree.
                o.insertRecord('trees',[i,t.getName(),t.getNewick(),t.getOrigin(),
                                        t.getScore()[0],t.getScore()[1],
                                        self.landscape.getNode(i)['explored']])
            else:
                print 'Warning: Tree %s has identical structure to %s.' % (
                    str(find),str(i))
        
        # Add the graph in its entirety as its adjacency list.
        adj_list = self.graph.edges_iter()
        for tupl in adj_list:
            source,target = [int(x) for x in tupl]
            o.insertRecord('graph',[source,target])
        
        # Add all of the locks.
        for l in self.landscape.getLocks():
            tree = l.topology.toTree()
            treeIndex = self.landscape.indexOf(tree)
            o.insertRecord('locks',[treeIndex,l.getBranchIndex()])
        
        # Close the file.
        o.close()

        # Return the path.
        return fpath

    def writeFile(self,path='.'):
        
        ''' Write the landscape serialized file to given path. '''

        return self._dump(path)

class landscapeParser(object):

    ''' Encapsulates the construction of a landscape
    object from a sqlite landscape file. '''

    def __init__(self,path):
        self.file = path
        self.metadata = None
        self.name = None
        self.database = None
        self.alignment = None
        self.treemap = {}
        self.trees = []
        self.explored = {}
        self.landscape = None

    def getName(self):
        
        ''' Acquire the name of the parsed landscape. '''
        
        if self.name is None:
            if len(self.metadata) > 0:
                for metadata in self.metadata:
                    if metadata[0] == 'name':
                        self.name = metadata[1]
            else: self.name = os.path.splitext(os.path.basename(self.file))[0]
        return self.name

    def _getLandscapeMetadata(self):
        
        ''' Check metadata in landscape if present. '''
        
        tables = self.database.getTables()
        if 'metadata' in tables: return self.database.getRecords('metadata')
        else: return []

    def _makeAlignment(self):
        
        ''' Construct a pseudo-FASTA string to create a new alignment from. '''
        
        pseudofasta = ''
        sequences = self.database.iterRecords('alignment')
        for sequence in sequences:
            _,name,seq = sequence
            if name != '' and seq != '':
                pseudofasta += '>%s\n%s\n' % (name,seq)
        if pseudofasta != '':
            self.alignment = alignment(str(pseudofasta))
        
    def _getTrees(self):
    
        floatIfNotNone = lambda d: float(d) if d != None else d
        intIfNotNone   = lambda d: int(d) if d != None else d
        for t in self.database.iterRecords('trees'):
            treeid,name,newick,orig,ml,pars,exp = t
            if newick != '':
                trobj = tree.tree(str(newick))
                trobj.origin = str(orig)
                trobj.name = name
                trobj.score = [floatIfNotNone(ml),intIfNotNone(pars)]
                self.treemap[int(treeid)] = trobj
                self.explored[trobj] = bool(exp)
        for t in sorted(self.treemap.keys()):
            self.trees.append(self.treemap[t])

    def _getGraph(self):
        
        getIDForTree = lambda d: self.landscape.findTreeTopologyByStructure(d.getStructure())
        for e in self.database.iterRecords('graph'):
            raw_s,raw_t = e
            source,target = getIDForTree(self.treemap[raw_s]),getIDForTree(self.treemap[raw_t])
            self.landscape.graph.add_edge(source,target)

    def _doLocks(self):

        for l in self.database.iterRecords('locks'):
            i,b = l
            self.landscape.lockBranchFoundInTreeByIndex(self.treemap[i],b)
        
    def parse(self):
        
        ''' Parse the file. '''

        # Get database object.
        fpath = self.file
        if not os.path.isfile(fpath): raise IOError('Landscape file not found.')
        o = SQLiteDatabase(fpath)
        self.database = o
        
        # Get metadata.
        self.metadata = self._getLandscapeMetadata()
        
        # Get the alignment.
        self._makeAlignment()
        
        # Get all of the trees.
        self._getTrees()
        
        # Construct the landscape.
        self.landscape = landscape(self.alignment,starting_tree=self.trees[0],
                                   operator=self.trees[0].getOrigin(),root=False)
        for t in self.trees:
            add = self.landscape.addTree(t)
            self.landscape.getVertex(add).setExplored(self.explored[t])
        
        # Ensure all edges are established.
        self._getGraph()
        
        # Apply locks.
        self._doLocks()

        return (self.landscape,self.getName())
