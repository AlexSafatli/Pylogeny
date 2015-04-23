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
import os, tree, sys

class landscapeWriter(object):

    ''' Encapsulate the writing of a landscape to a file format. '''

    def __init__(self, landscape, name):
        
        ''' Instantiates this writer.
        
        :param landscape: a landscape object
        :type landscape: a :class:`.landscape.landscape` object
        :param name: the name of this landscape
        :type name: a string
        
        '''
        
        # Fields
        self.landscape = landscape
        self.name      = name
        self.graph     = landscape.graph
        self.database  = None
        self.cleankey  = ''
        self.cleansuff = ''
        
        # Clean the input landscape name.
        self._cleanName()

    def _cleanName(self):

        if (type(self.name) != str):
            raise TypeError('Landscape name given not a string.')
        for x in self.name:
            # Replace all non-alphanumeric characters with a "_".
            if x.isalnum(): self.cleankey += x
            else: self.cleankey += '_'
        self.cleansuff = '%s.landscape' % (self.cleankey)

    def _schema(self):
        
        # Get the database object.
        dbobj = self.database
        
        # Create all content tables.
        dbobj.newTable('alignment',('seqid','integer'),('seqname','text'),
                       ('sequence','text'))
        dbobj.newTable('trees',('treeid','integer'),('name','text'),
                       ('newick','text'),('origin','text'),('ml','real'),
                       ('pars','real'),('explored','boolean'),
                       ('structure','text'))
        dbobj.newTable('graph',('source','integer'),('origin','integer'))
        dbobj.newTable('locks',('treeid','integer'),('branchid','integer'))
        
        # Create a table to store metadata.
        dbobj.newTable('metadata',('key','text'),('value','text'))

    def _dump(self,path='.'):

        # Open the file.
        fpath = os.path.join(path,self.cleansuff) # Generate file path.
        if os.path.isfile(fpath):
            os.unlink(fpath) # Remove file if already exists.
        o = SQLiteDatabase(fpath) # Open an SQLite database at that location.
        self.database = o # Assign this object to field.
        
        # Construct the schema for the landscape.
        self._schema()
        
        # Include metadata about the landscape. Required for parsing.
        o.insertRecord('metadata',['name',self.name])
        o.insertRecord('metadata',['version',VERSION])
        
        # Add the alignment as a set of sequence records labelled by taxa.
        index = 0
        if self.landscape.alignment != None:
            for s in self.landscape.alignment:
                o.insertRecord('alignment',[index,s.name,s.sequence])
                index += 1
        
        # Add all of the trees (most memory intensive elements in the DB).
        for i in self.landscape.iterNodes():
            t = self.landscape.getTree(i)
            s = t.getStructure()
            find = self.landscape.findTreeTopologyByStructure(s)
            if (find == None):
                raise AssertionError('Tree topology search failed (%s).' % (
                    str(i)))
            elif (find == i): # Verify unique identity of this tree.
                # Insert this tree as a record.
                name = t.getName()
                newi = t.getNewick()
                ori = t.getOrigin()
                scs = [t.getScore()[x] for x in xrange(0,len(t.getScore()))]
                o.insertRecord('trees',[i,name,newi,ori,scs[0],scs[1],
                                        self.landscape.getNode(i)['explored'],
                                        s])
            else:
                sys.stderr.write(
                    'Warning: Tree %s has identical structure to %s.\n' % (
                    str(find),str(i)))
        
        # Add the graph by copying its adjacency list.
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
        
        ''' Write the landscape serialized file to given path.
        
        :param path: a directory path, defaulting to the current one
        :type path: a string
        :return: the relative filepath to the written file
        
        '''

        return self._dump(path)

class landscapeParser(object):

    ''' Encapsulates the construction of a landscape
    object from a sqlite landscape file. '''

    def __init__(self,path):
        
        ''' Instantiate this parser. 
        
        :param path: the filepath to the landscape file
        :type path: a string
        
        '''
        
        self.file = path
        self.metadata = None
        self.name = None
        self.database = None
        self.alignment = None
        self.treemap = {}
        self.landscape = None

    def getName(self):
        
        ''' Acquire the name of the parsed landscape.
        
        :return: a string
        
        '''
        
        if self.name is None and self.metadata != None:
                for metadata in self.metadata:
                    if metadata[0] == 'name': self.name = metadata[1]
        elif self.name is None: self.name = os.path.splitext(
            os.path.basename(self.file))[0]
        return self.name

    def _getLandscapeMetadata(self):
        
        ''' Check metadata in landscape if present. '''
        
        tables = self.database.getTables()
        if 'metadata' in tables:
            return self.database.getRecords('metadata')
        else: raise IOError('Incompatible landscape DB format.')

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
            struct = None
            if len(t) == 7:
                treeid,name,newick,orig,ml,pars,exp = t
            else:
                treeid,name,newick,orig,ml,pars,exp,struct = t
            if newick != '':
                if (struct == None):
                    i = self.landscape.addTreeByNewick(str(
                        newick),score=False,check=False)
                else:
                    i = self.landscape.addTreeByNewick(str(
                        newick),score=False,check=False,struct=struct)
                trobj = self.landscape.getTree(i)
                trobj.setOrigin(str(orig))
                trobj.setName(name)
                trobj.setScore([floatIfNotNone(ml),intIfNotNone(pars)])
                self.treemap[int(treeid)] = trobj
                self.landscape.getVertex(i).setExplored(bool(exp))

    def _getGraph(self):
        
        getIDForTree = lambda d: self.landscape.findTreeTopologyByStructure(
            d.getStructure())
        for e in self.database.iterRecords('graph'):
            raw_s,raw_t = e
            source,target = getIDForTree(self.treemap[raw_s]),getIDForTree(
                self.treemap[raw_t])
            if source is None or target is None:
                raise IOError(
                    'Adjacency list record has source or target as None.')
            self.landscape.graph.add_edge(source,target)

    def _applyLocks(self):

        for l in self.database.iterRecords('locks'):
            self.landscape.lockBranchFoundInTreeByIndex(self.treemap[l[0]],l[1])
        
    def parse(self):
        
        ''' Parse the file.
        
        :return: a tuple of a :class:`.landscape.landscape` object and its name (a string)
        
        '''

        # Get database object.
        fpath = self.file
        if not os.path.isfile(fpath): raise IOError('Landscape file not found.')
        o = SQLiteDatabase(fpath)
        self.database = o
        
        # Get metadata.
        self.metadata = self._getLandscapeMetadata()
        
        # Get the alignment.
        self._makeAlignment()
        
        # Create an empty landscape structure.
        self.landscape = landscape(self.alignment,root=False)     
        
        # Add all of the trees.
        self._getTrees()
        
        # Set landscape operator.
        if len(self.landscape) > 0:
            index = self.landscape.getNodeNames()[-1]
            self.landscape.setOperator(self.landscape.getTree(
                index).getOrigin())
        
        # Ensure all edges are established.
        self._getGraph()
        
        # Apply locks.
        self._applyLocks()

        return (self.landscape,self.getName())
