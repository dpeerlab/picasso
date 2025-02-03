# Picasso

Picasso is a Python package for constructing phylogenetic trees from large-scale copy number alteration (CNA) data derived from single-cell sequencing. It implements an algorithm for inferring evolutionary relationships between cellular populations based on their CNA profiles.

## Installation

You can install Picasso using pip:
```bash
pip install git+https://github.com/dpeerlab/picasso.git
```

## Requirements

- Python 3.6+
- NumPy
- Pandas
- Seaborn
- Matplotlib

## Quick Start

```python
import picasso

# Load your CNA data
character_matrix = picasso.load_data()

# Create and fit the model
model = Picasso(character_matrix,
                min_depth=2,
                min_clone_size=5,
                assignment_confidence_threshold=0.8)
model.fit()

# Get the phylogeny and clone assignments
phylogeny = model.get_phylogeny()
clone_assignments = model.get_clone_assignments()
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

