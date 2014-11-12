// Date:   Oct 23 2014
// Author: Alex Safatli
// E-mail: safatli@cs.dal.ca

#include <Python.h>
#include <algorithm>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <locale>

using namespace std;

/* Phylogenetic Tree Implementation */

class Node {
  
  list<Node*> children;
  string parsimony_profile_data;
  string label;
  int len;
  
  public:

    // Constructors
    Node() {
      parsimony_profile_data = "";
      label = "";
      len = 0;
    };
  
    // Type Definitions
    typedef list<Node*> Children;
    typedef list<Node*>::iterator ChildIterator;
  
    // Function Definitions
    void setData(string str)          { parsimony_profile_data = str; };
    void setLabel(string str)         { label = str; };
    void addChild(Node *n)            { children.push_back(n); ++len; };
    string getData()                  { return parsimony_profile_data; };
    string getLabel()                 { return label; };
    Children *getChildren()           { return &children; };
    ChildIterator iterChildrenBegin() { return children.begin(); };
    ChildIterator iterChildrenEnd()   { return children.end();   };
    int numChildren()                 { return len; };
    bool isLeaf()                     { return (len == 0); };

};

class Tree {

  Node *root;

  public:

	  // Constructors/Destructor
	  Tree() { root = NULL; }
	  ~Tree() { deallocNode(root); }

	  // Type Definitions
    typedef Node* TreeNode;

    // Function Definitions
    TreeNode getRoot() { return root; };
    void setRoot(TreeNode n) { root = n; };
    TreeNode newNode() {
	    return new Node();
	  };
    vector<TreeNode> postOrder() {
	    vector<TreeNode> vec;
	    postOrder(root,&vec);
	    return vec;
	  };

  private:
  
    void deallocNode(TreeNode n) {
	    if (n != NULL) {
	      for (Node::ChildIterator it = n->iterChildrenBegin();
		    it != n->iterChildrenEnd(); ++it) deallocNode(*it);
        delete n;
      }
	  };
    void postOrder(TreeNode n, vector<TreeNode> *vec) {
	    if (n) {
        for (Node::ChildIterator it = n->iterChildrenBegin();
            it != n->iterChildrenEnd(); ++it) postOrder(*it,vec);
	      vec->push_back(n);
	    }
	  };

};

/* Newick String Parsing */

class PhylogenyParser {

  Tree *tree = NULL;
  string newick;
  int newicklength;

  public:

    // Constructor
    PhylogenyParser(Tree *tr, string s) { 
      tree = tr; newick = s;
      newicklength = newick.size();
    };

    // Function Definitions
    void parse() {
      parse(0,newick.size()-1,NULL);
    };

  private:

    int getBalancingBracket(int a) {
	    int opencount = 0, i = a;
	    while (i+1 < newicklength) {
	      if (newick[i+1] == '(') ++opencount;
	      else if (newick[i+1] == ')' && --opencount < 0) return i+1;
	      ++i;
	    } return -1;
	  };

    int readBranchLength(int a,float &f) {
	    int endPosOfFloat = a+1, i = a;
	    while (i+1 < newicklength) {
	      if (!isdigit(newick[i+1]) && newick[i+1] != '.') break;
	      else { ++endPosOfFloat; } ++i;
	    }
	    f = stof(newick.substr(a+1,endPosOfFloat-a-1));
	    return i;
	  };

	  int readLeafName(int a,string &name) {
	    int i = a;
	    while (i < newicklength) {
        if (newick[i] == ',' || newick[i] == ':' || newick[i] == ';' || 
            newick[i] == ')'|| i == newicklength) {
		      name.assign(newick,a,i-a);
	        return i-1;
	      } ++i;
      } return -1;
    };

    void parse(int i, int j, Tree::TreeNode cursor) {

	  // Parse a newick string into a topological newick structure given
	  // a top level node. Ported from Python code from in the Pylogeny
	  // python package. Input is substring [i,j] and cursor Node.

	  int q = -1, k = -1; float bl = 0.0;
	  while (i <= j) {

      // Create new node.
      Tree::TreeNode newnode = tree->newNode();

      // Add this new node as a child to the cursor node.
      if (cursor != NULL) cursor->addChild(newnode);
      else tree->setRoot(newnode); // Is root of tree.

	    // Two cases: (a) hit a subtree, or (b) hit a leaf.
	    if (newick[i] == '(') { // Subtree
		    k = getBalancingBracket(i);
		    if (newick[k+1] == ':') q = readBranchLength(k+1,bl);
		    else if (newick[k+1] != ',') {
	        string name = string();
	        q = readLeafName(k+1,name);
	        newnode->setLabel(name);
	        if (newick[q+1] == ':') q = readBranchLength(q+1,bl);
	      }
	      parse(i+1,k,newnode);
	      if (q != -1) k = q;
	    }
      else { // Leaf
		    string name = string();
		    k = readLeafName(i,name);
		    if (k+1<=j && newick[k+1] == ':') k = readBranchLength(k+1,bl);
		    newnode->setLabel(name);
	    }

      // Advance forward.
      i = k+1; // Next symbol.
      if (i<j) { while (i+1 <= j && newick[i] != ',') ++i; }
      ++i;

    }

  };

};

/* Intersection/Union Functions for String-Composed Sets */

string strIntersection(string a, string b) {

  string it;
  set_intersection(a.begin(),a.end(),b.begin(),b.end(),
                        back_inserter(it));
  return it;

}

string strUnion(string a, string b) {

  string it;
  set_union(a.begin(),a.end(),b.begin(),b.end(),
                back_inserter(it));
  return it;

}

string intersectChildrenDataOfNode(Node *n) {

  string interSet = string();
  for (Node::ChildIterator ch = n->iterChildrenBegin();
      ch != n->iterChildrenEnd(); ++ch) {
    string data = (*ch)->getData();
    if (interSet.size() == 0) interSet.assign(data);
    else interSet.assign(strIntersection(data,interSet));
  }
  return interSet;

}

string unionChildrenDataOfNode(Node *n) {

  string unionSet = string();
  for (Node::ChildIterator ch = n->iterChildrenBegin();
      ch != n->iterChildrenEnd(); ++ch) {
    string data = (*ch)->getData();
    if (unionSet.size() == 0) unionSet.assign(data);
    else unionSet.assign(strUnion(data,unionSet));
  }
  return unionSet;

}

/* Fitch Algorithm Implementation */

int cost(Tree *tree, int numPro, string *proArr, long *weiArr,
         map<string,int> *map) {

  // Calculate parsimony score. This works on multifurcating trees.
  int total = 0;
  vector<Tree::TreeNode> postorder = tree->postOrder();
  for (int i = 0; i < numPro; i++) {
    // For every profile...
    int local = 0;
    for (vector<Tree::TreeNode>::iterator it = postorder.begin();
         it != postorder.end(); it++) {
      if ((*it)->isLeaf()) {
        int taxInd = map->at((*it)->getLabel());
        (*it)->setData(proArr[i].substr(taxInd,1));
      }
      else {
        // Get the intersection of all children.
        string X = intersectChildrenDataOfNode(*it);
	      // Is intersection zero?
        if (X.size() == 0) {
          X = unionChildrenDataOfNode(*it);
          local += (*it)->numChildren() - 1; // Increase cost.
        } (*it)->setData(X);
      }
    }
    total += local*weiArr[i];
  }
  return total;

}

static PyObject * fitch_cost(PyObject *self, PyObject *args) {

  // Get input Newick string and lists of profiles & weights.
  char *newickstr; PyObject *proList, *weiList, *taxDict;
  if (!PyArg_ParseTuple(args,"sOOO",&newickstr,&proList,&weiList,&taxDict))
    return NULL;
  string newick(newickstr);

  // See how many elements comprise the list arguments.
  int proEles = PyList_Size(proList),
      weiEles = PyList_Size(weiList);
  if (proEles != weiEles) return NULL; // Do not match.

  // Construct and populate the tree.
  Tree tree;
  PhylogenyParser parser(&tree,newick);
  parser.parse();

  // Acquire profiles (first list argument).
  string *strArr = new string[proEles];
  for (int i = 0; i < proEles; i++) {
    char *proCurr; PyObject *curr;
    curr = PyList_GetItem(proList,i);
    if (PyString_Check(curr)) {
      proCurr = PyString_AsString(curr);
      strArr[i] = string(proCurr);
    } else return NULL;
  }

  // Acquire weights (second list argument).
  long *weiArr = new long[weiEles];
  for (int i = 0; i < proEles; i++) {
    PyObject *it = PyList_GetItem(weiList,i);
    if (PyInt_Check(it)) weiArr[i] = PyInt_AsLong(it);
    else return NULL;
  }

  // Acquire dictionary mapping taxa to vectors.
  map<string,int> taxMap;
  if (PyDict_Check(taxDict)) {
    PyObject *li = PyDict_Items(taxDict);
    int numEles = PyList_Size(li);
    for (int j = 0; j < numEles; j++) {
      PyObject *tu = PyList_GetItem(li,j);
      PyObject *key = PyTuple_GetItem(tu,0);
      PyObject *ind = PyTuple_GetItem(tu,1);
      string k = PyString_AsString(key);
      int v = (int)PyInt_AsLong(ind);
      taxMap.insert(pair<string,int>(k,v));
    }
  } else return NULL;

  // Get the fitch parsimony cost.
  int c = cost(&tree,proEles,strArr,weiArr,&taxMap);

  // Construct into a Python object and return.
  delete [] strArr;
  delete [] weiArr;
  return Py_BuildValue("i",c);

}

/* Python Extension Boilerplate */

static PyMethodDef modulemethods[] = {
  {"calculateCost",fitch_cost,METH_VARARGS,
  "Given a Newick string tree, calculate the parsimony cost."},
  {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initfitch(void) {
  (void) Py_InitModule("fitch",modulemethods);
}

int main(int argc, char *argv[]) {
  Py_SetProgramName(argv[0]);
  Py_Initialize();
  initfitch();
  return 0;
}
