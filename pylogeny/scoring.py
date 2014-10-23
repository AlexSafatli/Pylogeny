''' Phylogenetic tree scoring. '''

# Date:   Jan 29 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import pll, parsimony, fitch, p4
try:
    from model import DiscreteStateModel as State
    from pytbeaglehon.disc_state_cont_time_model import HKY85Model
    from pytbeaglehon.tree_scorer import create_toggle_partial_tree_scorer as Scorer
    from pytbeaglehon.tests.util import TreeForTesting as T
except: pass

# Scoring Functions

def beaglegetLogLikelihood(tree,alignment):
    
    ''' Acquire log-likelihood via C++ library
    BEAGLE via use of pybeaglethon wrapper library.
    Parameters: newick.tree object and alignment
    object. '''
    
    # Construct appropriate BEAGLE instance.
    if (alignment.getDim() == 4):
        stmod = alignment.getStateModel()
    else: return None
    model = HKY85Model(2.0,stmod.getStateFreqs())
    states = stmod.getAlignmentAsStateList()
    newi = tree.getSimpleNewick()[:-1]
    scorer = Scorer(model_list=[model],data=states,tree=T(newi))
    return scorer()

def getLogLikelihoodForTopology(topo,alignment):
    
    ''' Acquire log-likelihood via C library libpll.
    Parameters: rearrangement.topology object and alignment
    object. '''

    # Acquire a data model for libpll.
    try: p = pll.dataModel(topo,alignment)
    except: return None
    
    # Try and acquire the actual log-likelihood.
    try: sc = p.getLogLikelihood()
    except:
        p.close()
        return None

    p.close()
    return sc

def getLogLikelihood(tree,alignment):
    
    ''' Acquire log-likelihood via C library libpll.
    Parameters: newick.tree object and alignment
    object. '''

    return getLogLikelihoodForTopology(
        tree.toTopology(),alignment)

def getParsimonyFromProfiles(newick,profiles):
    
    ''' Acquire parsimony via a Python implementation.
    Parameters: Newick string and parsimony.profile_set object. '''
    
    topo = topology()
    topo.fromNewick(newick)
    return parsimony.fitch_cost(topo,profiles)

def getParsimony(newick,alignment):
    
    ''' Acquire parsimony via a Python implementation.
    Parameters: Newick string and alignment object. '''
     
    topo = topology()
    topo.fromNewick(newick)
    return parsimony.fitch(topo,alignment)

def getParsimonyFromProfilesForTopology(topology,profiles):

    ''' Acquire parsimony via a Python implementation. Parameters:
    rearrangement.topology object and alignment object. '''

    return parsimony.fitch(topology,profiles)

#def getParsimony(newick,alignment):
#    
#    ''' Acquire parsimony via a C++ implementation.
#    Parameters: newick string and alignment object. '''    
#    
#    profiles = parsimony.profile_set(alignment)
#    return getParsimonyFromProfiles(newick,profiles)

#def getParsimonyForTopology(topology,alignment):
#    
#    ''' Acquire parsimony via a C++ implementation.
#    Parameters: rearrangement.topology and 
#    alignment object. '''    
#    
#    profiles = parsimony.profile_set(alignment)
#    return getParsimonyFromProfiles(
#        topology.toNewick(),profiles)

#def getParsimonyFromProfiles(newick,profiles):
#    
#    ''' Acquire parsimony via a C++ implementation.
#    Parameters: newick string and parsimony.profile_set 
#    object. '''        
#    
#    prolist = [str(x) for x in profiles.profiles]
#    weilist = profiles.weights
#    taxdict = profiles.taxa
#    return fitch.calculateCost(newick,prolist,weilist,taxdict)

#def getParsimonyFromProfilesForTopology(topology,profiles):
#    
#    ''' Acquire parsimony via a C++ implementation.
#    Parameters: rearrangement.topology and parsimony.profile_set 
#    object. '''        
#    
#    prolist = [str(x) for x in profiles.profiles]
#    weilist = profiles.weights
#    taxdict = profiles.taxa
#    return fitch.calculateCost(topology.toNewick(),prolist,
#                               weilist,taxdict)
    
    
