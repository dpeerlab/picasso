API Reference
=============

This section contains the complete API documentation for PICASSO, automatically generated from the NumPy-style docstrings in the source code.

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   picasso.build_tree
   picasso.CloneTree
   picasso.utils
   picasso.itol_utils

Overview
--------

PICASSO (Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data Optimization) is designed to reconstruct tumor phylogenies from noisy, inferred copy number alteration (CNA) data derived from single-cell RNA sequencing.

Core Components
---------------

- :py:class:`picasso.build_tree.Picasso`: Main phylogenetic inference algorithm
- :py:class:`picasso.CloneTree.CloneTree`: Tree analysis and visualization tools  
- :py:mod:`picasso.utils`: Data preprocessing and loading utilities
- :py:mod:`picasso.itol_utils`: iTOL visualization export functions