# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2024-12-19

### Added
- Initial release of PICASSO (Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data Optimization)
- Core `Picasso` class for phylogenetic inference from noisy scRNA-seq-inferred CNA data
- `CloneTree` class for phylogenetic tree analysis and visualization
- Comprehensive utility functions for data preprocessing (`encode_cnvs_as_ternary`, `load_data`)
- iTOL export utilities for publication-quality visualizations
- Support for confidence-based termination criteria to handle noisy data
- Categorical mixture models for probabilistic clone assignment
- Comprehensive NumPy-style docstrings throughout
- Type hints for all public APIs
- Example dataset and Jupyter notebook tutorials

### Features
- **Noise-aware phylogenetic inference**: Specifically designed for scRNA-seq-inferred CNA data
- **Flexible termination criteria**: Both BIC and confidence-based stopping conditions
- **Comprehensive visualization**: Matplotlib-based plotting and iTOL export
- **Robust to noise**: Confidence thresholds and minimum clone sizes prevent over-fitting
- **Well-documented**: Extensive documentation with focus on noisy data handling
- **Type-safe**: Complete type annotation coverage

### Dependencies
- Python ≥ 3.10
- NumPy ≥ 1.24.0
- pandas ≥ 2.0.0
- pomegranate ≥ 1.0.4 (categorical mixture models)
- ete3 ≥ 3.1.3 (phylogenetic trees)
- matplotlib ≥ 3.6.0
- seaborn ≥ 0.12.0
- tqdm ≥ 4.64.0
- scipy ≥ 1.10.0

### Package Structure
```
picasso/
├── __init__.py          # Main package with comprehensive docstring
├── build_tree.py        # Core Picasso algorithm implementation
├── CloneTree.py         # Tree analysis and visualization
├── utils.py             # Data preprocessing utilities
├── itol_utils.py        # iTOL export functions
├── py.typed             # Type hint support marker
└── sample_data/         # Example datasets
```

### Installation
```bash
pip install picasso_phylo
# or
conda install -c conda-forge picasso_phylo
```

### Basic Usage
```python
from picasso import Picasso, CloneTree, load_data

# Load example data
cna_data = load_data()

# Reconstruct phylogeny
picasso = Picasso(cna_data, 
                 min_clone_size=10,
                 assignment_confidence_threshold=0.8)
picasso.fit()

# Analyze results
phylogeny = picasso.get_phylogeny()
assignments = picasso.get_clone_assignments()
clone_tree = CloneTree(phylogeny, assignments, cna_data)
clone_tree.plot_alterations()
```

## [Unreleased]

### Planned
- GPU acceleration for large datasets
- Additional termination criteria (chi-squared test)
- Enhanced visualization options
- Integration with scanpy ecosystem
- Benchmarking suite
- Advanced tutorial notebooks