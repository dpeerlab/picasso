Installation
============

PICASSO requires Python 3.10 or higher and has been tested on Linux, macOS, and Windows.

Quick Install
-------------

The easiest way to install PICASSO is via pip from PyPI:

.. code-block:: bash

   pip install picasso-phylo

The package is not (yet) available on conda-forge due to some dependency issues. To use it in a conda or mamba environment, please install via pip inside your environment:

.. code-block:: bash

   conda create -n picasso_env python=3.10
   conda activate picasso_env
   pip install picasso-phylo

System Requirements
-------------------

**Python Version**
   - Python 3.10, 3.11, or 3.12
   - Older Python versions are not supported due to type hint requirements

**Operating Systems**
   - Linux (Ubuntu 18.04+, CentOS 7+)
   - macOS (10.14+) 
   - Windows 10/11

**Memory Requirements**
   - Minimum: 4 GB RAM
   - Recommended: 8+ GB RAM for datasets with >10,000 cells
   - Large datasets (>50,000 cells) may require more RAM

**Dependencies**
   PICASSO automatically installs the following required packages:

   - ``numpy >= 1.21.0`` - Numerical computing
   - ``pandas >= 1.3.0`` - Data manipulation and analysis
   - ``ete3 >= 3.1.2`` - Phylogenetic tree handling
   - ``pomegranate >= 1.0.0`` - Mixture model implementation
   - ``scipy >= 1.7.0`` - Scientific computing utilities
   - ``tqdm >= 4.62.0`` - Progress bars
   - ``matplotlib >= 3.5.0`` - Basic plotting (optional for visualization)
   - ``seaborn >= 0.11.0`` - Statistical plotting (optional for visualization)

Installation Methods
--------------------

Standard Installation (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For most users, the standard pip installation is sufficient:

.. code-block:: bash

   pip install picasso-phylo

This installs PICASSO and all required dependencies.

Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~

If you want to contribute to PICASSO or need the latest development features:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/dpeerlab/picasso.git
   cd picasso

   # Install in development mode
   pip install -e .

   # Or install with development dependencies
   pip install -e ".[dev]"

Virtual Environment Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend using a virtual environment (conda/mamba or venv) to avoid dependency conflicts.

.. code-block:: bash

   # Create conda environment
   conda create -n picasso-env python=3.11
   conda activate picasso-env
   
   # Use pip within conda environment
   pip install picasso-phylo


Verification
------------

Test your installation by running:

.. code-block:: python

   import picasso
   print(f"PICASSO version: {picasso.__version__}")
   
   # Load example data
   data = picasso.load_data()
   print(f"Example data shape: {data.shape}")
   
   # Quick test run
   model = picasso.Picasso(data, min_clone_size=10)
   print("Installation successful!")

Expected output:

.. code-block:: text

   PICASSO version: 1.0.0
   Example data shape: (10000, 116)
   Installation successful!

Troubleshooting
---------------

**Common Issues**

*Import Error: No module named 'picasso'*
   - Ensure you installed ``picasso-phylo``, not ``picasso`` (different package)
   - Check you're using the correct Python environment
   - Try: ``pip install --upgrade picasso-phylo``

*Pomegranate Installation Issues*
   - Some systems may require: ``pip install cython numpy`` before installing PICASSO
   - For M1/M2 Macs: ``conda install pomegranate`` may work better than pip

*Performance Issues with Large Datasets*
   - Reduce dataset size or increase system memory. Filter features with low variance before analysis.
   - Use larger ``min_clone_size`` parameter to reduce computational complexity
   - Consider running on a high-memory system or cluster
   - Use appropriate termination criteria to avoid over-fitting

**Getting Help**

If you encounter issues not covered here:

1. Check the `GitHub Issues <https://github.com/dpeerlab/picasso/issues>`_ page
2. Search existing issues or create a new one
3. Include your Python version, operating system, and error message
4. Provide a minimal example that reproduces the problem

Next Steps
----------

Once installation is complete, proceed to the :doc:`quickstart` guide to learn basic usage, or explore the :doc:`examples` for more detailed tutorials.