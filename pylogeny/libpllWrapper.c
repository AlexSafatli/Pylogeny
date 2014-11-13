/* libpllWrapper.c
   Date: Jan 31, 2014
   Author: Alex E. Safatli
----------------------------------
   Uncomprehensive bindings with libpll - phylogenetic functions
   intended for tree scoring and rearrangement. Designed for libpll 1.0.0
*/

#include <pll/pll.h>
#include <Python.h>
#include <stdlib.h>

/* Problem Structure */

typedef struct {
   pllAlignmentData   *alignment;
   pllInstanceAttr      *attribs;
   pllInstance         *instance;
   pllNewickTree           *tree;
   partitionList     *partitions;
   pllRearrangeList *arrangelist;
} problem_t;

/* Function Prototypes */

static PyObject * phylo_rearr_to_pylist(problem_t * probl);
static PyObject * phylo_init(PyObject *self, PyObject *args);
static PyObject * phylo_del(PyObject *self, PyObject *args);
static PyObject * phylo_getml(PyObject *self, PyObject *args);
static PyObject * phylo_getspr_withindist(PyObject *self, PyObject *args);
static PyObject * phylo_getnni_withindist(PyObject *self, PyObject *args);
static PyObject * phylo_getall_withindist(PyObject *self, PyObject *args);
PyMODINIT_FUNC initlibpllWrapper(void);
int main(int argc, char *argv[]);

/* Python Object Wrapping */

static PyObject * phylo_rearr_to_pylist(problem_t * probl) {

   pllRearrangeList *rlist; PyObject *list; int i;
   rlist = probl->arrangelist;
   list  = PyList_New(rlist->entries);

   /* Iterate over items in rearrangement list, commit the
      rearrangement, get the Newick string, store other information
      into list entry, rollback. */

   for (i = 0; i < rlist->entries; i++) {

      // Set up data structure(s).
      PyObject * tupl;
      char * rtype = (rlist->rearr[i].rearrangeType
                      == PLL_REARRANGE_SPR) ? "SPR" : "NNI";
      char * rnewick;
      float  rlikeli = rlist->rearr[i].likelihood;

      // Commit the move.
      pllRearrangeCommit(probl->instance,probl->partitions,
                         &(rlist->rearr[i]),PLL_TRUE);
      Tree2String(probl->instance->tree_string,probl->instance,
                  probl->partitions,probl->instance->start->back,
                  PLL_TRUE,PLL_TRUE,0,0,0,PLL_SUMMARIZE_LH,0,0);

      // Get the newick string.
      rnewick = probl->instance->tree_string;

      // Rollback the move.
      pllRearrangeRollback(probl->instance,probl->partitions);

      // Flatten to Python tuple and put into list.
      tupl = PyTuple_Pack(3,Py_BuildValue("s",rtype),
                            Py_BuildValue("f",rlikeli),
                            Py_BuildValue("s",rnewick));
      PyList_SetItem(list,i,tupl);

   }

   pllEvaluateGeneric(probl->instance,probl->partitions,
                      probl->instance->start,PLL_TRUE,PLL_FALSE);
   return list;

}

/* Component Functions */

static PyObject * phylo_init(PyObject *self, PyObject *args) {

   // Declarations
   const char  *alignf,
               *newick,
               *partf;
   problem_t   *problem;

   // Parse arguments.
   if (!PyArg_ParseTuple(args,"sss",&alignf,&newick,&partf)) {
      PyErr_SetString(PyExc_IOError,"Could not acquire correct input.");
      return NULL;
   }

   // Allocate memory appropriately from heap.
   problem = (problem_t*)malloc(sizeof(problem_t));
   if (!problem) {
      PyErr_SetString(PyExc_ReferenceError,"Memory allocation failure.");
      return NULL;
   }

   // Create attribute set for new problem.
   pllInstanceAttr *attr  = problem->attribs
      = (pllInstanceAttr*)malloc(sizeof(pllInstanceAttr));
   if (!problem->attribs) {
      PyErr_SetString(PyExc_ReferenceError,"Memory allocation failure.");
      return NULL;
   }
      
   // Set attributes.
   attr->rateHetModel     =  PLL_GAMMA;
   attr->fastScaling      =  PLL_FALSE;
   attr->saveMemory       =  PLL_FALSE;
   attr->useRecom         =  PLL_FALSE;
   attr->randomNumberSeed = 0xFACEFACE;
   attr->numberOfThreads  =          8;

   // Create new problem instance.
   problem->instance = pllCreateInstance(attr);

   // Parse files, Newick string.
   problem->alignment = pllParseAlignmentFile(PLL_FORMAT_PHYLIP,alignf);
   problem->tree      = pllNewickParseString(newick);
   pllQueue *partinfo = pllPartitionParse(partf);
   if (!problem->alignment || !problem->tree) {
      if (!problem->alignment)
         PyErr_SetString(PyExc_IOError,"Alignment file could not be parsed.");
      else PyErr_SetString(PyExc_IOError,"Tree was not able to be parsed.");
      return NULL;
   }

   // Commit partitions; build structure up.
   pllAlignmentData *ali = problem->alignment;
   pllInstance     *inst = problem->instance;
   problem->partitions   = pllPartitionsCommit(partinfo,ali);
   pllTreeInitTopologyNewick(inst,problem->tree,PLL_TRUE);
   if (!pllLoadAlignment(inst,ali,problem->partitions,PLL_DEEP_COPY)) {
      PyErr_SetString(PyExc_IOError,"Could not find correct correspondances.");
      return NULL;
   }

   // Start model.
   pllInitModel(inst,problem->partitions,ali);
   pllQueuePartitionsDestroy(&partinfo);

   /* Setup rearrangement list. Maximum number of rearrangements
      will be the maximum number of SPR neighbors: 4(n-3)(n-2) + (n-1). */
   int n = inst->mxtips, max = 4 * (n-3) * (n-2) + (n-1);
   problem->arrangelist = pllCreateRearrangeList(max);

   // Tell libpll to optimize length of new branch created by moves.
   inst->thoroughInsertion = PLL_TRUE;

   // Return pointer to problem.
   return PyCObject_FromVoidPtr((void*)problem,NULL);

}

static PyObject * phylo_del(PyObject *self, PyObject *args) {

   // Get pointer.
   PyObject *p; problem_t *pp;
   if (!PyArg_ParseTuple(args,"O",&p) || !PyCObject_Check(p)) {
      PyErr_SetString(PyExc_ReferenceError,"Could not parse pointer.");
      return NULL;
   } pp = (problem_t*)PyCObject_AsVoidPtr(p);

   // Deallocate everything. Not incl. the universe.
   pllDestroyRearrangeList(&(pp->arrangelist));
   pllAlignmentDataDestroy(pp->alignment);
   pllNewickParseDestroy(&(pp->tree));
   pllPartitionsDestroy(pp->instance,&(pp->partitions));
   pllDestroyInstance(pp->instance);
   free(pp->attribs); free(pp);

   // Return None.
   return Py_BuildValue("");

}

static PyObject * phylo_getml(PyObject *self, PyObject *args) {

   // Get pointer.
   PyObject *p; problem_t *pp;
   if (!PyArg_ParseTuple(args,"O",&p) || !PyCObject_Check(p)) {
      PyErr_SetString(PyExc_ReferenceError,"Could not parse pointer.");
      return NULL;
   } pp = (problem_t*)PyCObject_AsVoidPtr(p);

   // Get log-likelihood.
   pllTreeEvaluate(pp->instance,pp->partitions,64);

   // Return it.
   if (pp->instance->likelihood)
      return Py_BuildValue("f",pp->instance->likelihood);
   else {
      PyErr_SetString(PyExc_IOError,"Unable to evaluate tree for likelihood.");    
      return NULL;
   }

}

static PyObject * phylo_getspr_withindist(PyObject *self, PyObject *args) {

   // Get pointer, radius to go to.
   PyObject *p; int max_dist; problem_t *pp;
   if (!PyArg_ParseTuple(args,"Oi",&p,&max_dist) || !PyCObject_Check(p)) {
      PyErr_SetString(PyExc_ReferenceError,"Could not parse pointer, integer.");
      return NULL;
   } pp = (problem_t*)PyCObject_AsVoidPtr(p);

   int n = pp->instance->mxtips, i;

   for (i = 1; i <= n+1; i++) {
      pllRearrangeSearch(pp->instance,pp->partitions,PLL_REARRANGE_SPR,
                         pp->instance->nodep[i],1,max_dist,
                         pp->arrangelist);
   }

   return phylo_rearr_to_pylist(pp);

}

static PyObject * phylo_getnni_withindist(PyObject *self, PyObject *args) {

   // Get pointer, radius to go to.
   PyObject *p; int max_dist; problem_t *pp;
   if (!PyArg_ParseTuple(args,"Oi",&p,&max_dist) || !PyCObject_Check(p)) {
      PyErr_SetString(PyExc_ReferenceError,"Could not parse pointer, integer.");
      return NULL;
   } pp = (problem_t*)PyCObject_AsVoidPtr(p);

   int n = pp->instance->mxtips, i;

   for (i = 1; i <= n+1; i++) {
      pllRearrangeSearch(pp->instance,pp->partitions,PLL_REARRANGE_NNI,
                         pp->instance->nodep[i],1,max_dist,
                         pp->arrangelist);
   }

   return phylo_rearr_to_pylist(pp);

}

static PyObject * phylo_getall_withindist(PyObject *self, PyObject *args) {

   // Get pointer, radius to go to.
   PyObject *p; int max_dist; problem_t *pp;
   if (!PyArg_ParseTuple(args,"Oi",&p,&max_dist) || !PyCObject_Check(p)) {
      PyErr_SetString(PyExc_ReferenceError,"Could not parse pointer, integer.");
      return NULL;
   } pp = (problem_t*)PyCObject_AsVoidPtr(p);

   int n = pp->instance->mxtips, i;

   /*
   pllRearrangeSearch(pp->instance,pp->partitions,
                      PLL_REARRANGE_SPR,
                      pp->instance->nodep[n+1],1,max_dist,
                      pp->arrangelist);

   pllRearrangeSearch(pp->instance,pp->partitions,
                      PLL_REARRANGE_NNI,
                      pp->instance->nodep[n+1],1,max_dist,
                      pp->arrangelist); */
   
   for (i = 1; i <= n+1; i++) {
      pllRearrangeSearch(pp->instance,pp->partitions,PLL_REARRANGE_NNI,
                         pp->instance->nodep[i],1,max_dist,
                         pp->arrangelist);
      pllRearrangeSearch(pp->instance,pp->partitions,PLL_REARRANGE_SPR,
                         pp->instance->nodep[i],1,max_dist,
                         pp->arrangelist);
   }

   return phylo_rearr_to_pylist(pp);

}

/* Python Extension Boilerplate */

static PyMethodDef modulemethods[] = {
   {"new", phylo_init, METH_VARARGS,
   "Initialize a problem with an alignment, Newick tree, and partition model."},
   {"destroy", phylo_del, METH_VARARGS,
   "Destroy an input problem and deallocate all resources."},
   {"getLogLikelihood", phylo_getml, METH_VARARGS,
   "Score a tree and acquire its log-likelihood for a given input problem."},
   {"getSPRMovesInDistance", phylo_getspr_withindist, METH_VARARGS,
   "Get all SPR moves within a maximum distance from a leaf of current given tree instance."},
   {"getNNIMovesInDistance", phylo_getnni_withindist, METH_VARARGS,
   "Get all NNI moves within a maximum distance from a leaf of current given tree instance."},
   {"getAllMovesInDistance", phylo_getall_withindist, METH_VARARGS,
   "Get all moves within a maximum distance from a leaf of current given tree instance."},
   {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initlibpllWrapper(void) {
   (void) Py_InitModule("libpllWrapper",modulemethods);
}

int main(int argc, char *argv[]) {
   Py_SetProgramName(argv[0]);
   Py_Initialize();
   initlibpllWrapper();
   return 0;
}
