PICASSO: Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data Optimization
======================================================================================================

.. image:: https://img.shields.io/pypi/v/picasso-phylo.svg
   :target: https://pypi.org/project/picasso-phylo/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/pyversions/picasso-phylo.svg
   :target: https://pypi.org/project/picasso-phylo/
   :alt: Python versions

.. image:: https://readthedocs.org/projects/picasso-phylo/badge/?version=latest
   :target: https://picasso-phylo.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://img.shields.io/github/license/yourusername/picasso.svg
   :target: https://github.com/yourusername/picasso/blob/main/LICENSE
   :alt: License

**PICASSO** is a Python package for reconstructing tumor phylogenies from noisy, inferred copy number alteration (CNA) data derived from single-cell RNA sequencing (scRNA-seq). Unlike methods designed for direct single-cell DNA sequencing data, PICASSO is specifically optimized to handle the uncertainty and noise inherent in CNAs inferred from gene expression profiles.

Key Features
------------

**Noise-Aware Phylogeny Reconstruction**
   - Handles uncertainty in scRNA-seq-inferred CNA data
   - Probabilistic assignment with confidence thresholds
   - Robust to technical artifacts and dropout events

**Flexible Tree Building**
   - Iterative binary splitting with categorical mixture models
   - Multiple termination criteria (BIC, confidence-based, chi-squared)
   - Customizable depth and clone size constraints

**Comprehensive Analysis**
   - Clone aggregation and modal profile generation
   - Evolutionary change inference along tree branches
   - Integration with iTOL for publication-ready visualizations

**Designed for Single-Cell Data**
   - Optimized for the specific challenges of scRNA-seq CNA inference
   - Handles variable clone sizes and imbalanced datasets
   - Prevents over-fitting to noise patterns

Quick Start
-----------

Installation
~~~~~~~~~~~~

Install PICASSO from PyPI:

.. code-block:: bash

   pip install picasso-phylo

Or from conda-forge:

.. code-block:: bash

   conda install -c conda-forge picasso-phylo

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from picasso import Picasso, load_data

   # Load example CNA data
   cna_data = load_data()

   # Initialize PICASSO with noise-appropriate parameters
   picasso = Picasso(cna_data, 
                    min_clone_size=10,  # Larger for noisy data
                    assignment_confidence_threshold=0.8)

   # Reconstruct phylogeny
   picasso.fit()

   # Get results
   phylogeny = picasso.get_phylogeny()
   clone_assignments = picasso.get_clone_assignments()

   # Analyze results
   from picasso import CloneTree
   tree_analyzer = CloneTree(phylogeny, clone_assignments, cna_data)
   tree_analyzer.plot_alterations()

Algorithm Overview
------------------

PICASSO uses an iterative binary splitting approach:

1. **Initialization**: All cells start in a single root clone
2. **Iterative Splitting**: At each depth level:
   
   - Fit Categorical Mixture Models with k=1 and k=2 components for each clone
   - Evaluate splitting criteria (BIC or assignment confidence)
   - Split clones that meet the criteria into two daughter clones

3. **Termination**: Stop when no clones meet splitting criteria or constraints are reached
4. **Tree Construction**: Build phylogenetic tree from the clone hierarchy

The algorithm is specifically designed to handle:

- **Noise and artifacts** in scRNA-seq-inferred CNAs
- **Uncertainty** in copy number state assignments  
- **Variable clone sizes** and imbalanced data
- **Over-fitting** to noise patterns through confidence-based termination

Citation
--------

If you use PICASSO in your research, please cite our preprint while our manuscript is under review:

.. code-block:: bibtex

   @article{picasso2025,
   title={Transcriptomic plasticity is a hallmark of metastatic pancreatic cancer},
   author={Jim{'e}nez-S{'a}nchez, Alejandro and Persad, Sitara and Hayashi, Akimasa and Umeda, Shigeaki and Sharma, Roshan and Xie, Yubin and Mehta, Arnav and Park, Wungki and Masilionis, Ignas and Chu, Tinyi and Zhu, Feiyang and Hong, Jungeui and Chaligne, Ronan and O'Reilly, Eileen M. and Mazutis, Linas and Nawy, Tal and Pe'er, Itsik and Iacobuzio-Donahue, Christine A. and Pe'er, Dana},
   journal={bioRxiv},
   year={2025},
   month={February},
   day={28},
   doi={10.1101/2025.02.28.640922},
   url={https://www.biorxiv.org/content/10.1101/2025.02.28.640922v1},
   note={Preprint}
   }
   
Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   examples

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api

.. toctree::
   :maxdepth: 1
   :caption: Development
   
   changelog

Support
-------

- **Documentation**: https://picasso-phylo.readthedocs.io/
- **Source Code**: https://github.com/yourusername/picasso
- **Issue Tracker**: https://github.com/yourusername/picasso/issues
- **PyPI Package**: https://pypi.org/project/picasso-phylo/

License
-------

PICASSO is released under the MIT License. See the `LICENSE <https://github.com/yourusername/picasso/blob/main/LICENSE>`_ file for details.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`