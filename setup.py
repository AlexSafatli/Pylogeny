''' An installation script for Pylogeny. '''

# Date:   Oct 16 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from setuptools import setup, Extension as extension
import os

# Metadata

VERSION = '0.2.3.1'
DESCRIP = 'A code framework for phylogenetic tree reconstruction, rearrangement, scoring, and for the manipulation, heuristic search, and analysis of the phylogenetic tree combinatorial space.'
LONG    = 'A Python library and code framework for phylogenetic tree reconstruction, rearrangement, scoring, and for the manipulation and analysis of the phylogenetic tree combinatorial space. Also possesses features to execute popular heuristic programs such as FastTree and RAxML to acquire approximate ML trees.'
URL     = 'http://www.github.com/AlexSafatli/Pylogeny'
AUTHOR  = 'Alex Safatli'
EMAIL   = 'safatli@cs.dal.ca'
DEPNDS  = ['numpy','networkx','pandas']
LINKS   = ['http://p4-phylogenetics.googlecode.com/archive/4491de464e68fdb49c7a11e06737cd34a98143ec.tar.gz']
PKGDATA = {'pylogeny':['fitch.cpp','pylibpll.c']}
#FITCHCC = os.path.join('pylogeny','fitch.cpp')
PLLC    = os.path.join('pylogeny','pylibpll.c')

# Compilation for C/C++ Extensions (Fitch, Pylibpll)

pllExtension   = extension('pylibpll',sources=[PLLC],include_dirs=['/usr/local/include'],libraries=['pll-generic'],library_dirs=['/usr/local/lib'])
#fitchExtension = extension('fitch',sources=[FITCHCC],include_dirs=)

# Setup

setup(name='pylogeny',version=VERSION,description=DESCRIP,long_description=LONG,url=URL,author=AUTHOR,author_email=EMAIL,license='MIT',packages=['pylogeny'],package_data=PKGDATA,ext_modules=[pllExtension],install_requires=DEPNDS,dependency_links=LINKS,zip_safe=False)
