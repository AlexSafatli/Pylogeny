''' Pylogeny is a Python library and code framework for phylogenetic tree reconstruction and scoring in the larger frame of a combinatorial space.

Allows one to perform the following tasks: (1) Generate and manage phylogenetic tree landscapes. (2) Build and rearrange phylogenetic trees using preset operators such as NNI, SPR, and TBR. (3) Score phylogenetic trees by Log-likelihood and Parsimony. 

Dependencies: Numpy, NetworkX, Pandas, P4 Phylogenetic Library, libpll. Suggested: FastTree, RAxML, PytBEAGLEhon. '''

import os as _

thisdir = _.path.split(_.path.realpath(__file__))[0]
itlist = _.listdir(thisdir)
__all__ = [_.path.split(x)[-1].strip('.py') for x in itlist if x.endswith('.py') and not x.endswith('__init__.py')]
