"""
PICASSO: Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data Optimization.

PICASSO is a computational method for reconstructing tumor phylogenies from noisy,
inferred copy number alteration (CNA) data derived from single-cell RNA sequencing
(scRNA-seq). Unlike methods designed for direct scDNA-seq data, PICASSO specifically
handles scalability and the uncertainty and noise inherent in CNA profiles inferred from gene
expression data.

Key Features
------------
- **Noise-aware phylogenetic inference**: Uses probabilistic models to handle
  uncertainty in scRNA-seq-inferred CNAs
- **Confidence-based termination**: Prevents over-fitting to noise through
  assignment confidence thresholds
- **Flexible model selection**: Supports both BIC and confidence-based criteria
  for clone splitting decisions
- **Comprehensive visualization**: Integrated plotting and iTOL export capabilities
- **Scalable implementation**: Handles datasets with hundreds to thousands of cells

Main Components
---------------
Picasso : class
    Core phylogenetic inference algorithm with categorical mixture modeling
CloneTree : class
    Phylogenetic tree analysis and visualization framework
utils : module
    Data preprocessing utilities including ternary encoding
itol : module
    Export functions for iTOL-compatible visualizations

Quick Start
-----------
>>> import pandas as pd
>>> from picasso import Picasso, CloneTree, load_data
>>>
>>> # Load example noisy CNA data
>>> cna_data = load_data()
>>>
>>> # Reconstruct phylogeny with noise-appropriate parameters
>>> picasso = Picasso(cna_data,
...                  min_clone_size=10,  # Larger for noisy data
...                  assignment_confidence_threshold=0.8,
...                  terminate_by='probability')
>>> picasso.fit()
>>>
>>> # Extract results
>>> phylogeny = picasso.get_phylogeny()
>>> assignments = picasso.get_clone_assignments()
>>>
>>> # Create integrated analysis object
>>> clone_tree = CloneTree(phylogeny, assignments, cna_data)
>>> clone_tree.plot_alterations(save_as='cna_heatmap.pdf')

Typical Workflow
----------------
1. **Data preparation**: Load CNA matrix (cells × genomic features)
2. **Parameter selection**: Choose parameters appropriate for noise level
3. **Phylogenetic inference**: Run PICASSO algorithm with `fit()`
4. **Result extraction**: Get phylogeny and clone assignments
5. **Visualization**: Use CloneTree or iTOL export for publication figures
6. **Analysis**: Examine clonal evolution patterns and CNA profiles

Parameters for Noisy Data
-------------------------
For very noisy scRNA-seq-inferred CNAs, consider:
- `min_clone_size=20-50`: Larger clone sizes reduce spurious clusters
- `assignment_confidence_threshold=0.85-0.95`: Higher confidence requirements
- `assignment_confidence_proportion=0.9-0.95`: More stringent cell proportions
- `max_depth=6-10`: Limit tree depth to prevent over-fitting
- `bic_penalty_strength=1.2-2.0`: Stronger penalty for model complexity

Installation Requirements
-------------------------
- Python ≥ 3.10
- pomegranate ≥ 1.0 (categorical mixture models)
- pandas, numpy (data handling)
- ete3 (phylogenetic trees)
- matplotlib, seaborn (visualization)
- tqdm (progress tracking)

Citation
--------
If you use PICASSO in your research, please cite our paper.

See Also
--------
- iTOL web interface: https://itol.embl.de/
- ETE Toolkit documentation: http://etetoolkit.org/
- Pomegranate documentation: https://pomegranate.readthedocs.io/
"""

# Core classes and functions
from .build_tree import Picasso
from .CloneTree import CloneTree
from .utils import encode_cnvs_as_ternary, load_data

# Import itol_utils as submodule for organized access
from . import itol_utils as itol

# Version information
__version__ = "0.1.0"
__author__ = "Sitara Persad, Alejandro Jimenez-Sanchez"
__email__ = "sitara.persad@columbia.edu"

# Define public API
__all__ = [
    # Core classes
    "Picasso",
    "CloneTree",
    # Utility functions
    "encode_cnvs_as_ternary",
    "load_data",
    # Submodules
    "itol",
    # Package metadata
    "__version__",
    "__author__",
    "__email__",
]
