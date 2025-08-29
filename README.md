# PICASSO: Phylogenetic Inference of Copy number Alterations in Single-cell Sequencing data Optimization

[![PyPI version](https://badge.fury.io/py/picasso-phylo.svg)](https://badge.fury.io/py/picasso-phylo)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/picasso_phylo.svg)](https://anaconda.org/conda-forge/picasso_phylo)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/github/license/dpeerlab/picasso)](LICENSE)

PICASSO is a computational method for reconstructing tumor phylogenies from noisy, inferred copy number alteration (CNA) data derived from single-cell RNA sequencing (scRNA-seq). Unlike methods designed for direct scDNA-seq data, PICASSO specifically handles the uncertainty and noise inherent in CNA profiles inferred from gene expression data.

## Key Features

- **Noise-aware phylogenetic inference**: Uses probabilistic models to handle uncertainty in scRNA-seq-inferred CNAs
- **Confidence-based termination**: Prevents over-fitting to noise through assignment confidence thresholds  
- **Comprehensive visualization**: Integrated plotting and iTOL export capabilities
- **Scalable implementation**: Handles datasets with hundreds to thousands of cells
- **Well-documented**: Extensive documentation with focus on noisy data handling

## Installation

### PyPI (recommended)
```bash
pip install picasso_phylo
```

### Conda
The package is not (yet) available on conda-forge due to some dependency issues. To use it in a conda or mamba environment, please install via pip inside your environment:
```bash
conda create -n picasso_env python=3.10
conda activate picasso_env
pip install picasso-phylo
```

### Development Installation
```bash
git clone https://github.com/dpeerlab/picasso
cd picasso
pip install -e ".[dev]"
```

## Requirements

- **Python**: ≥ 3.10
- **Core dependencies**: numpy, pandas, pomegranate, ete3, matplotlib, seaborn, tqdm, scipy
- **Optional**: jupyter (notebooks), pyqt5 (advanced visualization)

## Quick Start

```python
from picasso import Picasso, CloneTree, load_data

# Load example CNA data
cna_data = load_data()

# Initialize PICASSO with noise-appropriate parameters
picasso = Picasso(cna_data,
                 min_clone_size=10,  # Larger for noisy data
                 assignment_confidence_threshold=0.8,
                 terminate_by='probability')

# Reconstruct phylogeny
picasso.fit()

# Extract results
phylogeny = picasso.get_phylogeny()
assignments = picasso.get_clone_assignments()

# Create integrated analysis object
clone_tree = CloneTree(phylogeny, assignments, cna_data)
clone_tree.plot_alterations(save_as='cna_heatmap.pdf')
```

### For Very Noisy scRNA-seq Data

```python
# Use stricter parameters for very noisy data
picasso_strict = Picasso(cna_data,
                        min_clone_size=50,
                        max_depth=8,  # Limit depth
                        assignment_confidence_threshold=0.9,
                        assignment_confidence_proportion=0.95,
                        bic_penalty_strength=1.5)
picasso_strict.fit()
```

## Features

### Data Processing
- Load and process copy number alteration (CNA) data
- Encode CNVs as ternary values for more meaningful similarity measures
- Feature selection to remove non-informative regions

### Tree Construction
- Construct phylogenetic trees using the PICASSO algorithm
- Flexible tree manipulation and rooting options
- Support for both clone-level and sample-level phylogenies

### Visualization
- Basic tree visualization
- Clone size plotting
- Alteration plotting
- Integration with iTOL for advanced visualization
- Support for:
  - Heatmaps
  - Colorstrips
  - Stacked bar charts

## Advanced Usage

### Tree Construction and Manipulation

```python
from picasso import CloneTree

# Create and manipulate the clone tree
tree = CloneTree(phylogeny, clone_assignments, filtered_matrix, clone_aggregation='mode')
outgroup = tree.get_most_ancestral_clone()
tree.root_tree(outgroup)

# Get different tree representations
clone_tree = tree.get_clone_phylogeny()
cell_tree = tree.get_sample_phylogeny()
```

### iTOL Visualization

```python
# Generate heatmap of copy number changes
heatmap_annot = picasso.itol.dataframe_to_itol_heatmap(character_matrix)
with open('heatmap_annotation.txt', 'w') as f:
    f.write(heatmap_annot)

# Generate colorstrip annotation
colorstrip_annot = picasso.itol.dataframe_to_itol_colorstrip(
    data_series,
    color_map,
    dataset_label='Label'
)

# Generate stacked bar visualization
stackedbar_annot = picasso.itol.dataframe_to_itol_stackedbar(
    proportions_df,
    color_map,
    dataset_label='Label'
)
```

## API Reference

### Picasso Class Parameters

- `min_depth`: Minimum depth of the phylogenetic tree
- `max_depth`: Maximum depth of the tree (None for unlimited)
- `min_clone_size`: Minimum number of samples in a clone
- `terminate_by`: Criterion for terminating tree growth
- `assignment_confidence_threshold`: Confidence threshold for sample assignment
- `assignment_confidence_proportion`: Required proportion of samples meeting confidence threshold
- `bic_penalty_strength`: Strength of BIC penalty term. Higher values (>1.0) encourage simpler models, useful for noisy data to prevent over-fitting.

## Visualization

For detailed visualization, we recommend using the [iTOL website/application](https://itol.embl.de/), which accepts newick strings as input and allows for detailed customization of tree visualization. Picasso provides convenience functions for generating iTOL annotation files to visualize metadata on the tree.

## Support

If you encounter any problems, please open an issue along with a detailed description.

## License

This project is licensed under the MIT License:

```
MIT License

Copyright (c) 2024 [Pe'er Lab]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Citation

If you use Picasso in your research, please cite our paper. 

