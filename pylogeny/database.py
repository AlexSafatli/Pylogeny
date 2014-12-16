''' Connect, access, + manipulate external tree data from a remote SQL server or from a sqlite file. '''

# Date:   Nov 7 2013
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import MySQLdb as mysql, sqlite3 as sqllite
from landscape import landscape, choice
from tree import tree
from abc import ABCMeta as abstractclass, abstractmethod

class DatabaseLandscape(landscape):
    
    ''' Abstract the landscape to one comprising a landscape. '''
    
    __metaclass__ = abstractclass

    def _acquireNumberOfTaxa(self):
        
        topo = self.getTree(0).toTopology()
        self.leaves = len(topo.getLeaves())
        
    @abstractmethod
    def _fetchTreeFromDatabase(self,i): pass
    
    @abstractmethod
    def _fetchRearrangementsFromDatabase(self,i): pass

    def getNode(self,i):
        
        self._fetchTreeFromDatabase(self,i)
        return super(DatabaseLandscape,self).getNode(self,i)

class SQLExhaustiveLandscape(DatabaseLandscape):
    
    def __init__(self,dbobj,aliname):
        
        # Set up fields.
        super(SQLLandscape,self).__init__(None,None,False)
        self.database = dbobj
        self.aliname  = aliname
        self.queryID  = lambda d: dbobj.filterRecords(
            'newick','ID=%d' % (d)) # Alias lookup function.
        self.queryScr = lambda d: dbobj.filterRecords(
            aliname,'TID=%d' % (d))
        self.tables   = dbobj.getTables()
        if not aliname in self.tables:
            raise IOError('Could not find alignment table <%s>.' % (aliname))
        self._fetchTreeFromDatabase(0) # Get root tree.
        self._acquireNumberOfTaxa()
    
    def _fetchTreeFromDatabase(self,i):
        
        # Overrides lookup in graph object first. Looks in database.
        if i in self.graph.node: return
        entries = self.queryID(i)
        if len(entries) == 0: return   
        entry = entries[0]
        
        # Create tree object.
        t = tree(entry[1])
        t.setName(i)
        self.addTree(t) 
        
        # Get relevant score(s).
        entry = self.queryScr(i)[0]
        ml,p  = entry[2:4]
        t.score = (float(ml),float(p))
    
    def _fetchRearrangementsFromDatabase(self,i):
        moves = self.database.filterRecords('operations','TID=%d'%(i))
        return [(int(x),int(y),_,int(z)) for x,y,_,z in moves]
        
    def _applyRearrangement(self,move):
        e,v,o = move[1],move[-1],move[2]
        self.graph.add_edge(e,v)
        self.getTree(v).origin = o
        
    def exploreRandomTree(self,i):
        node  = self.getNode(i)
        if (node['explored']): return None
        moves = self._fetchRearrangementsFromDatabase(i)
        while (len(moves) > 0):
            move = choice(moves)
            moves.remove(move)
            j = move[-1]
            if j in self.graph.node and self.graph.has_edge(i,j):
                continue
            self._fetchTreeFromDatabase(j)
            self._applyRearrangement(move)
            self.getEdge(i,j)['weight'] = self.defaultWeight
            return j
        node['explored'] = True
    
    def getDatabaseNode(self,i):
        self._fetchTreeFromDatabase(i)
        return super(SQLLandscape,self).getNode(i)
    
    def exploreTree(self,i):
        node  = self.getNode(i)
        if (node['explored']): return None
        moves = self._fetchRearrangementsFromDatabase(i)
        while (len(moves) > 0):
            move = choice(moves)
            moves.remove(move)
            j = move[-1]
            if j in self.graph.node and self.graph.has_edge(i,j):
                continue
            self._fetchTreeFromDatabase(j)
            self._applyRearrangement(move)
            self.getEdge(i,j)['weight'] = self.defaultWeight
        node['explored'] = True     
        return True

class SQLiteLandscape(landscape):
    
    ''' Allow random access of all landscape data from an sqlite file found on the 
    hard disk. '''
    
    def __init__(self,dbobj):
        
        pass
        

class database(object):

    ''' Allow interfacing with a SQL/sqlite database. '''
    
    __metaclass__ = abstractclass
    cursor = None # Cursor for the database object.
    
    @abstractmethod
    def getTables(self): pass
    
    @abstractmethod
    def getColumns(self,table): pass  
    
    def isEmpty(self):
        ''' Determine if the database is empty. '''
        return (len(self.getTables()) == 0)
    def getHeaders(self,table):
        ''' Get only header names for a given table's columns. '''
        return [x[0] for x in self.getColumns(table)]
    def getRecordsColumn(self,table,col):
        ''' Get all data for a single colmun from records for a table. '''
        self.query("""SELECT %s FROM %s""" % (col,table))
        return self.cursor.fetchall()
    def getRecords(self,table):
        ''' Get all records from a given table in the database. '''
        self.query("""SELECT * FROM %s""" % (table))
        return self.cursor.fetchall()
    def iterRecords(self,table):
        ''' Get a record, one at a time, from a table in the database. '''
        self.query("""SELECT * FROM %s""" % (table))
        nextitem = self.cursor.fetchone()
        while nextitem != None:
            yield nextitem
            nextitem = self.cursor.fetchone()
    def filterRecords(self,table,condn):
        ''' Get all records from a given table following a condition. '''
        self.query("""SELECT * FROM %s WHERE %s""" % (table,condn))
        return self.cursor.fetchall()
    def getRecordsAsDict(self,table):
        ''' Acquires records using getRecords() and then leverages
        access using a dictionary data structure. '''
        d = {}
        items = self.getRecords(table)
        if len(items) == 0: return {}
        headers = self.getHeaders(table)
        if len(items[-1]) != len(headers): return None
        for header in headers: d[header] = []
        for item in items:
            for i in xrange(len(headers)):
                header = headers[i]
                value  = item[i]
                d[header].append(value)
        return d
    def newTable(self,tablename,*args):
        self.query("""CREATE TABLE %s (%s)""" % (
            tablename,', '.join(['%s %s' % (x,y) for x,y in args])))
    def insertRecords(self,tablename,items):
        self.querymany('''INSERT INTO %s VALUES (%s)''' % (
            tablename,', '.join(['?' for x in items[0]])),items)
    def insertRecord(self,tablename,record):
        self.insertRecords(tablename,[record])

    @abstractmethod
    def query(self,q): pass
    
    @abstractmethod
    def querymany(self,q,i): pass

    @abstractmethod
    def close(self): pass

class SQLDatabase(database):
    
    ''' Database object to allow reading from a MySQL database. '''
    
    def __init__(self,host,user,pw,db):
        self.hostname = host
        self.username = user
        self.password = pw
        self.database = db
        self.connect()
    
    def connect(self):
        self.socket   = mysql.connect(host=self.hostname,
                                      user=self.username,
                                      passwd=self.password,
                                      db=self.database)
        self.cursor   = self.socket.cursor()        
        
    def getTables(self):
        self.query("""SHOW TABLES""")
        return [x[0] for x in self.cursor.fetchall()]        
        
    def getColumns(self,table):
        ''' Return column information for a given table. '''
        self.query("""SHOW COLUMNS FROM %s""" % (table))
        return self.cursor.fetchall()          
        
    def query(self,q):
        # Execute a MySQL query.
        try:
            self.cursor.execute(q)
        except (AttributeError,mysql.OperationalError):
            # In case of server timeout, reconnect.
            self.connect()
            self.query(q)        

    def querymany(self,q,i):
        try:
            self.cursor.executemany(q,i)
        except (AttributeError,mysql.OperationalError):
            self.connect()
            self.querymany(q,i)

    def close(self):
        # Close the connection.
        self.socket.close()
        

class SQLiteDatabase(database):
    
    def __init__(self,filepath):
        self.filepath = filepath
        self.socket   = sqllite.connect(filepath)
        self.cursor   = self.socket.cursor()

    def getColumns(self,table):
        ''' Return column information for a given table. '''
        self.query("""PRAGMA table_info(%s)""" % (table))
        return self.cursor.fetchall()

    def getTables(self):
        self.query("""SELECT name FROM sqlite_master WHERE type='table';""")
        return [x[0] for x in self.cursor.fetchall()]   

    def query(self,q):
        # Execute an sqlite query.
        self.cursor.execute(q)

    def querymany(self,q,i):
        self.cursor.executemany(q,i)

    def close(self):
        # Close the connection.
        self.socket.commit()
        self.socket.close()    

        
