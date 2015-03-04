Pylogeny
========

[![Build Status](https://travis-ci.org/AlexSafatli/Pylogeny.svg?branch=master)](https://travis-ci.org/AlexSafatli/Pylogeny)

A software library and code framework, written in the Python programming language for Python 2.6+, for phylogenetic tree reconstruction, rearrangement, scoring, and for the manipulation, heuristic search, and analysis of the phylogenetic tree combinatorial space. Functionality also exists in the framework to execute popular heuristic programs such as FastTree and RAxML to acquire approximate ML trees.

The following tasks are capable of being performed with this library:

  - Generate and maintain phylogenetic tree landscapes.
  - Construct and analyse heuristic methods to search these spaces.
  - Build and rearrange phylogenetic trees using preset operators (such as NNI, SPR, and TBR).
  - Score phylogenetic trees by Maximum Likelihood (calculated as log-likelihood) and Parsimony.
  - Build confidence sets of trees using the widely known [CONSEL](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/ "CONSEL") application.

[[ Note that C and C++ extensions are not built to work with Python 3.3+. You are advised to use a version 2 of the Python interpreter to run this package. ]]

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

You can print the Newick string for the global optimum or tree with the maximum likelihood by calling the following function on that object.

    globalMax = ls.getGlobalOptimum()
    print ls.getTree(globalMax)

Installation
-------------

You must first install non-Python dependencies, which currently stand as only [libpll](http://libpll.org). Acquire the appropriate binary or install, from source, version 1.0.2 with SSE3 support.

Once you have acquired and installed all of the necessary non-Python dependencies, you can install this software automatically using `pip` or `easy_install` with the command

    pip install pylogeny

or the command

    easy_install pylogeny

Documentation
-------------

The API is currently being output to a Github page located [here](http://AlexSafatli.github.io/Pylogeny "Pylogeny API").

Dependencies
-------------

 * [NumPy](http://www.numpy.org/)
 * [NetworkX](https://networkx.github.io/)
 * [Pandas](http://pandas.pydata.org/)
 * MySQLdb for Python
 * [P4](https://code.google.com/p/p4-phylogenetics/) Phylogenetic Library
 * [libpll](http://libpll.org) Phylogenetic Likelihood Library (currently compiled with latest [version 1.0.2 with SSE3](http://libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz))

Works With
-------------

 * [FastTree](http://www.microbesonline.org/fasttree/)
 * [RAxML](http://sco.h-its.org/exelixis/software.html)
 * [CONSEL](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/)
 * [PytBEAGLEhon](https://github.com/mtholder/pytbeaglehon)

Contributing
-------------

To contribute to this project, feel free to make a pull request and it will be reviewed by the code maintainers.