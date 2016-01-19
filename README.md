Pylogeny
========

A software library and code framework, written in the Python programming language for Python 2, for phylogenetic tree reconstruction, rearrangement, scoring, and for the manipulation, heuristic search, and analysis of the phylogenetic tree combinatorial space. Scoring of trees in this library is accomplished by bindings to the [libpll](http://libpll.org) phylogenetic C library. Functionality also exists in the framework to execute popular heuristic programs such as FastTree and RAxML to acquire approximate ML trees.

The following tasks are capable of being performed with this library:

  - Generate and maintain phylogenetic tree landscapes.
  - Construct and analyse heuristic methods to search these spaces.
  - Build and rearrange phylogenetic trees using preset operators (such as NNI, SPR, and TBR).
  - Score phylogenetic trees by Maximum Likelihood (calculated as log-likelihood) and Parsimony.
  - Build confidence sets of trees using the widely known [CONSEL](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/ "CONSEL") application.

Code Example
-------------

You can create a landscape for a given sequence alignment, add a tree to the landscape corresponding to the one acquired from FastTree, and then perform a hill-climbing search on that landscape on the basis of parsimony with the below code.

    from pylogeny.alignment import alignment
    from pylogeny.landscape import landscape
    from pylogeny.heuristic import parsimonyGreedy

    ali = alignment('yourAlignment.fasta')
    ls  = landscape(ali,starting_tree=ali.getApproxMLTree())
    heu = parsimonyGreedy(ls,ls.getRootNode())
    heu.explore()     

This performs an exploration but the heuristic does not return anything. In order to acquire an idea of what this fitness landscape now looks like, you can start making queries to the landscape.

    for tree in ls.iterTrees():
        print tree # Print all the Newick strings in the landscape.
    globalMax = ls.getGlobalOptimum() # Get the name for tree with best score.
    print ls.getTree(globalMax) # Print this tree (see its Newick string).

We can also see what neighboring trees have been explored from the first tree we started the heuristic at.

    neighbors = ls.getNeighborsFor(ls.getRoot())
    for neighbor in neighbors: # Neighbors are indices.
        print ls.getTree(neighbor) # Print the Newick string for that neighbor.

Installation
-------------

Installation requires access to a UNIX-like system or terminal. Furthermore, **basic build tools**, **python development header include files**, and [**MySQL library bindings**](https://www.mysql.com/) are required to install this software. In Ubuntu, this is done with the command

     sudo apt-get install python-dev build-essential libmysqlclient-dev

If you do not use a Debian or Ubuntu-derived Linux distribution, search for instructions on acquiring these for your platform.

Before continuing, the non-Python library dependency [libpll](http://libpll.org) must be installed. Acquire the appropriate binary or build, from source, version 1.0.2 with SSE3 support. A **convenient shell script** is located in the root directory of this repository that will perform this installation. You can run this without any download by peforming the command:

    wget https://raw.githubusercontent.com/AlexSafatli/Pylogeny/master/install-pll.sh -O - | sh

Once you have acquired and installed all of the necessary non-Python dependencies, you can install this software automatically using `pip` or `easy_install` with command

    pip install pylogeny

or

    easy_install pylogeny

respectively.

Documentation
-------------

Generated documentation is found [here](http://AlexSafatli.github.io/Pylogeny "Pylogeny API"). Tutorials and additional code examples are present in the [wiki](https://github.com/AlexSafatli/Pylogeny/wiki) on GitHub for this project.

Dependencies
-------------

 * [NetworkX](https://networkx.github.io/) version >= 1.9.1
 * [Pandas](http://pandas.pydata.org/) version >= 0.15.2
 * MySQLdb for Python version >= 1.2.5
 * [P4](https://code.google.com/p/p4-phylogenetics/) Phylogenetic Library version >= 0.93
 * [libpll](http://libpll.org) Phylogenetic Likelihood Library (currently compiled with latest [version 1.0.2 with SSE3](http://libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz))

Works With
-------------

 * [FastTree](http://www.microbesonline.org/fasttree/)
 * [RAxML](http://sco.h-its.org/exelixis/software.html)
 * [CONSEL](http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/consel/)
 * [PytBEAGLEhon](https://github.com/mtholder/pytbeaglehon)

Citing
-------------

To cite this library, refer to the paper *Pylogeny: an open-source Python framework for phylogenetic tree reconstruction and search space heuristics* that has been published to describe its purpose [here](https://peerj.com/articles/cs-9/) at [PeerJ](https://peerj.com/).

Contributing
-------------

To contribute to this project, feel free to make a pull request and it will be reviewed by the code maintainers.
