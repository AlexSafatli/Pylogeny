#include <Gulo/NewickFormat>
#include <Python.h>
#include <algorithm>
#include <vector>
#include <string>
#include <map>

using namespace Gulo;

/* Phylogenetic Tree Implementation & Gulo Glue */

class Node {
  // Binary tree node ONLY.
  std::string profile_data;
  std::string label;
  Node *left;
  Node *right;
  public:
    Node() { profile_data = ""; label = ""; left = NULL; right = NULL; };
    void setData(std::string str) { profile_data = str; };
    void setLabel(std::string str) { label = str; };
    void addChild(Node *n) { if (left == NULL) left = n;
                             else if (right == NULL) right = n; };
    std::string getData() { return profile_data; };
    std::string getLabel() { return label; };
    Node *getLeft() { return left; };
    Node *getRight() { return right; };
    bool isLeaf() { return (left == NULL && right == NULL); };
};

class Tree {
  Node *root;
  std::vector<Node*> li;
  public:
    typedef Node* TreeNode;
    typedef std::string TreeEdge;
    Tree() { root = NULL; };
    ~Tree() { freeNode(root); };
    TreeNode getRoot() { return root; };
    void setRoot(TreeNode n) { root = n; };
    TreeNode addNode() { TreeNode n = new Node(); li.push_back(n); return n; };
    std::vector<TreeNode> postOrder();
  private:
    void freeNode(TreeNode leaf);
    void postOrder(TreeNode n, std::vector<TreeNode> *vec);
};

void Tree::freeNode(Tree::TreeNode leaf) {
  if (leaf != NULL) {
    freeNode(leaf->getLeft());
    freeNode(leaf->getRight());
    delete leaf;
  }
}

std::vector<Tree::TreeNode> Tree::postOrder() {
  std::vector<Tree::TreeNode> v;
  postOrder(root,&v);
  return v;
}

void Tree::postOrder(Tree::TreeNode n, std::vector<Tree::TreeNode> *vec) {
  if (n) {
    postOrder(n->getLeft(),vec);
    postOrder(n->getRight(),vec);
    vec->push_back(n);
  }
}

namespace Gulo {

  template <>
  struct PhylogenyParserMapping<NewickFormat,Tree> {

    typedef typename Tree::TreeNode NodeIdentifier;
    typedef typename Tree::TreeEdge EdgeIdentifier;

    NodeIdentifier addNode(Tree &t, const std::vector<std::string> &comments) {
      return t.addNode();
    }

    EdgeIdentifier addEdge(Tree &t, NodeIdentifier from, NodeIdentifier to) {
      from->addChild(to);
      return from->getLabel();
    }

    EdgeIdentifier trailingEdge(Tree&) {
      return "";
    }

    void setLength(Tree &t, EdgeIdentifier edge, double len, const std::vector<std::string> &comments) {
      ;
    }

    void setLabel(Tree &t, NodeIdentifier node, std::string &&label, const std::vector<std::string> &comments) {
      node->setLabel(label);
    }

    void setHybridNodeType(Tree &t, NodeIdentifier n, const std::string &label, const std::vector<std::string> &comments) {
      ;
    }

    void setRoot(Tree &t, NodeIdentifier no) {
      t.setRoot(no);
    }

    void setUnrooted(Tree &t) {
      ;
    }

    void abort(Tree &t) {
      ;
    }

  };

}

/* Intersection/Union Functions for String-Composed Sets */

std::string st_intersection(std::string a, std::string b) {

  std::string it;
  std::set_intersection(a.begin(),a.end(),b.begin(),b.end(),
                        std::back_inserter(it));
  return it;

}

std::string st_union(std::string a, std::string b) {

  std::string it;
  std::set_union(a.begin(),a.end(),b.begin(),b.end(),
                std::back_inserter(it));
  return it;

}

/* Fitch Algorithm Implementation */

int cost(Tree *tree, int numPro, std::string *proArr, long *weiArr,
         std::map<std::string,int> *map) {

  // Calculate parsimony score.
  int total = 0;
  std::vector<Node*> postorder = tree->postOrder();
  for (int i = 0; i < numPro; i++) {
    // For every profile...
    int local = 0;
    for (std::vector<Node*>::iterator it = postorder.begin();
         it != postorder.end(); it++) {
      if ((*it)->isLeaf()) {
        int taxInd = map->at((*it)->getLabel());
        (*it)->setData(proArr[i].substr(taxInd,1));
      }
      else {
        // Get both sets.
        Node *a, *b;
        a = (*it)->getLeft();
        b = (*it)->getRight();
        std::string sa,sb;
        sa = a->getData();
        sb = b->getData();
        // Get intersection.
        std::string X = st_intersection(sa,sb);
        if (X.length() == 0) {
          local++; X = st_union(sa,sb);
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
  std::string newick(newickstr);

  // See how many elements comprise the list arguments.
  int proEles = PyList_Size(proList),
      weiEles = PyList_Size(weiList);
  if (proEles != weiEles) return NULL; // Do not match.

  // Make a tree object using Boost library.
  Tree tree;
  PhylogenyParser<NewickFormat,Tree> parser;
  parser.parse(newick,tree);

  // Acquire profiles (first list argument).
  std::string strArr[proEles];
  for (int i = 0; i < proEles; i++) {
    char *proCurr; PyObject *curr;
    curr = PyList_GetItem(proList,i);
    if (PyString_Check(curr)) {
      proCurr = PyString_AsString(curr);
      strArr[i] = std::string(proCurr);
    } else return NULL;
  }

  // Acquire weights (second list argument).
  long weiArr[weiEles];
  for (int i = 0; i < proEles; i++) {
    PyObject *it = PyList_GetItem(weiList,i);
    if (PyInt_Check(it)) weiArr[i] = PyInt_AsLong(it);
    else return NULL;
  }

  // Acquire dictionary mapping taxa to vectors.
  std::map<std::string,int> taxMap;
  if (PyDict_Check(taxDict)) {
    PyObject *li = PyDict_Items(taxDict);
    int numEles = PyList_Size(li);
    for (int j = 0; j < numEles; j++) {
      PyObject *tu = PyList_GetItem(li,j);
      PyObject *key = PyTuple_GetItem(tu,0);
      PyObject *ind = PyTuple_GetItem(tu,1);
      std::string k = PyString_AsString(key);
      int v = (int)PyInt_AsLong(ind);
      taxMap.insert(std::pair<std::string,int>(k,v));
    }
  } else return NULL;

  // Get the fitch parsimony cost.
  int c = cost(&tree,proEles,strArr,weiArr,&taxMap);

  // Construct into a Python object and return.
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
