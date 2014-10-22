Pylogeny
========

A Python library and code framework for phylogenetic tree reconstruction, rearrangement, scoring, and for the manipulation, heuristic search, and analysis of the phylogenetic tree combinatorial space. Also possesses features to execute popular heuristic programs such as FastTree and RAxML to acquire approximate ML trees.

Allows one to perform the following tasks:

  - Generate and maintain phylogenetic tree landscapes.
  - Build and rearrange phylogenetic trees using preset operators such as NNI, SPR, and TBR.
  - Score phylogenetic trees by Log-likelihood and Parsimony.
  - Build confidence sets of trees using CONSEL.

Code Example
-------------

You can create a landscape for a given sequence alignment, add a tree to the landscape corresponding to the one acquired from FastTree, and then perform a hill-climbing search on that landscape on the basis of parsimony with the below code.

    from pylogeny.alignment import alignment
    from pylogeny.landscape import landscape
    from pylogeny.heuristic import parsimonyGreedy

    ali = alignment('yourAlignment.fasta')
    ls  = landscape(ali,starting_tree=ali.getApproxMLTree())
    heu = parsimonyGreedy(ls,ls.getRootTree())
    heu.explore()     

Dependencies
-------------

 * Numpy
 * NetworkX
 * Pandas
 * P4 Phylogenetic Library
 * libpll Phylogenetic Likelihood Library

Suggested
-------------

 * FastTree
 * RAxML
 * CONSEL
 * PytBEAGLEhon
