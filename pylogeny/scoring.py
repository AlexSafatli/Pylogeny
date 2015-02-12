''' Functions for phylogenetic tree goodness-of-fit scoring. '''

# Date:   Jan 29 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import pll, parsimony, fitch, p4, tree, rearrangement, alignment
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
    Currently uses HKY85 model. 
    
    :param tree: A tree object.
    :type tree: :class: `tree.tree`
    :param alignment: An alignment object.
    :type alignment: :class: `alignment.alignment`
    :returns: A floating point value.
    
    '''
    
    # Construct appropriate BEAGLE instance.
    if (alignment.getDim() == 4):
        stmod = alignment.getStateModel()
    else: return None
    model = HKY85Model(2.0,stmod.getStateFreqs())
    states = stmod.getAlignmentAsStateList()
    newi = tree.getSimpleNewick()[:-1]
    scorer = Scorer(model_list=[model],data=states,tree=T(newi))
    return scorer()

def getLogLikelihood(tree,alignment,updateBranchLengths=True):
    
    ''' Acquire log-likelihood via C library libpll. 

    :param tree: A tree object.
    :type tree: :class: `tree.tree`
    :param alignment: An alignment object.
    :type alignment: :class: `alignment.phylipFriendlyAlignment`
    :param updateBranchLengths: Whether or not to update the branch lengths in the provided tree with optimized ones.
    :returns: A floating point value.    
    
    '''

    topo = tree.toTopology()
    
    # Acquire a data model for libpll.
    try: p = pll.dataModel(topo,alignment)
    except: return None

    # Try and acquire the actual log-likelihood.
    try: sc = p.getLogLikelihood()
    except:
        p.close()
        return None
    
    # Update the tree's Newick string based on new branch lengths.
    if updateBranchLengths:
        try: tree.updateNewick(p.getNewickString().strip(),reroot=True)
        except ValueError, e: pass
    
    # Close instance and return.
    p.close()
    return sc    

def _pygetParsimonyFromProfiles(topology,profiles):
    
    ''' Acquire parsimony via a Python implementation.
    Parameters: rearrangement.topology object and 
    parsimony.profile_set object. Deprecated function. '''
    
    return parsimony.fitch_cost(topology,profiles)

def _pygetParsimony(topology,alignment):
    
    ''' Acquire parsimony via a Python implementation.
    Parameters: rearrangement.topology object and 
    alignment object. Deprecated function. '''
     
    return parsimony.fitch(topology,alignment)

def getParsimony(newick,alignment):
    
    ''' Acquire parsimony via a C++ implementation.
    
    :param newick: A New Hampshire (Newick) tree string.
    :param alignment: An alignment object.
    :type alignment: :class: `alignment.alignment`
    :returns: An integer value.
    
    '''    
    
    profiles = parsimony.profile_set(alignment)
    return getParsimonyFromProfiles(newick,profiles)

def getParsimonyForTopology(topo,alignment):
    
    ''' Acquire parsimony via a C++ implementation.
    
    :param topo: A topology object.
    :type topo: :class: `rearrangement.topology`
    :param alignment: An alignment object.
    :type alignment: :class: `alignment.alignment`
    :returns: An integer value.
    
    '''    
    
    profiles = parsimony.profile_set(alignment)
    return getParsimonyFromProfiles(
        topology.toNewick(),profiles)

def getParsimonyFromProfiles(newick,profiles):
    
    ''' Acquire parsimony via a C++ implementation.
    
    :param newick: A New Hampshire (Newick) tree string.
    :param profiles: A set of profiles corresponding to an alignment.
    :type profiles: :class: `parsimony.profile_set`
    :returns: An integer value.    

    '''        
    
    prolist = [str(x) for x in profiles.profiles]
    weilist = profiles.weights
    taxdict = profiles.taxa
    return fitch.calculateCost(newick,prolist,weilist,taxdict)

def getParsimonyFromProfilesForTopology(topology,profiles):
    
    ''' Acquire parsimony via a C++ implementation.
    
    :param topo: A topology object.
    :type topo: :class: `rearrangement.topology`
    :param profiles: A set of profiles corresponding to an alignment.
    :type profiles: :class: `parsimony.profile_set`
    :returns: An integer value.
    
    '''        
    
    prolist = [str(x) for x in profiles.profiles]
    weilist = profiles.weights
    taxdict = profiles.taxa
    return fitch.calculateCost(topology.toNewick(),prolist,
                               weilist,taxdict)
    
    
