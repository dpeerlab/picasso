Detailed Examples
=================

This section provides comprehensive examples for different types of analyses with PICASSO. Each example demonstrates best practices for specific scenarios you might encounter with single-cell CNA data.

.. contents:: Examples in this guide
   :local:
   :depth: 2

Example 1: Basic Phylogeny Reconstruction
------------------------------------------

This example demonstrates the standard PICASSO workflow using the built-in example dataset.

**Setup and Data Loading**

.. code-block:: python

   import picasso
   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   
   # Load example CNA data
   character_matrix = picasso.load_data()
   print(f'Dataset: {character_matrix.shape[0]} cells × {character_matrix.shape[1]} features')
   
   # Examine the data structure
   print("Data range:", character_matrix.min().min(), "to", character_matrix.max().max())
   print("First few rows:")
   print(character_matrix.head())

**Data Visualization**

.. code-block:: python

   # Visualize the CNA data as a heatmap
   plt.figure(figsize=(12, 8))
   sns.heatmap(character_matrix.iloc[:100],
               cmap='coolwarm', center=0,
               cbar_kws={'label': 'Copy Number'})
   plt.title('Copy Number Alterations (first 100 cells)')
   plt.xlabel('Genomic Features')
   plt.ylabel('Cells')
   plt.tight_layout()
   plt.show()

**Data Encoding for Complex CNAs**

.. code-block:: python

   # Encode complex CNAs as ternary values for better similarity handling
   character_matrix = picasso.encode_cnvs_as_ternary(character_matrix)
   
   # Visualize the encoded data
   plt.figure(figsize=(12, 8))
   sns.heatmap(character_matrix.iloc[:100],
               cmap='coolwarm', center=0,
               cbar_kws={'label': 'Copy Number'})
   plt.title('Copy Number Alterations (first 100 cells)')
   plt.xlabel('Genomic Features')
   plt.ylabel('Cells')
   plt.tight_layout()
   plt.show()

**Basic Phylogenetic Reconstruction**

.. code-block:: python

   # Initialize PICASSO with standard parameters
   model = picasso.Picasso(
       character_matrix,
       min_clone_size=5,
       assignment_confidence_threshold=0.8,
       assignment_confidence_proportion=0.9
   )
   
   # Fit the model
   print("Reconstructing phylogeny...")
   model.fit()
   
   # Extract results
   phylogeny = model.get_phylogeny()
   clone_assignments = model.get_clone_assignments()
   
   print(f"Reconstructed phylogeny with {len(phylogeny.get_leaves())} terminal clones")

**Visualize Clone Size Distribution**

.. code-block:: python

   # Plot the distribution of clone sizes
   plt.figure(figsize=(10, 6))
   sns.ecdfplot(clone_assignments['clone_id'].value_counts())
   plt.title('Distribution of Clone Sizes')
   plt.xlabel('Clone Size')
   plt.ylabel('Count')
   plt.show()

**Alternative: BIC-based Termination**

.. code-block:: python

   # We can also use BIC-based termination; for small datasets, 
   # it may terminate with less resolved clones
   model = picasso.Picasso(
       character_matrix,
       min_clone_size=5,
       terminate_by='BIC'
   )
   
   # Fit the model
   print("Reconstructing phylogeny...")
   model.fit()
   
   # Extract results
   phylogeny = model.get_phylogeny()
   clone_assignments = model.get_clone_assignments()
   
   print(f"Reconstructed phylogeny with {len(phylogeny.get_leaves())} terminal clones")
   print(f"Clone size distribution:")
   print(clone_assignments['clone_id'].value_counts().head())

**Tree Analysis & Downstream Visualization**

.. code-block:: python

   # Create CloneTree for advanced analysis
   tree = picasso.CloneTree(phylogeny, clone_assignments, character_matrix)
   
   # Root the tree at the most ancestral clone
   outgroup = tree.get_most_ancestral_clone()
   tree.root_tree(outgroup)
   print(f"Tree rooted at clone: {outgroup}")
   
   # Generate visualizations showing the clones and their groupings (not phylogenetic structure)
   tree.plot_clone_sizes()
   tree.plot_alterations()
   
   # Get clone phylogeny as Newick string for external tools
   clone_tree = tree.get_clone_phylogeny()
   print("Newick format (first 100 characters):")
   print(clone_tree.write()[:100] + "...")

Example 2: Filtering Very Noisy scRNA-seq Data
-----------------------------------------------

This example shows how to handle very noisy CNA data typically obtained from scRNA-seq inference.

**Data Preparation**

.. code-block:: python

   # Load data
   character_matrix = picasso.load_data()
   
   # Encode complex CNAs as ternary values for better similarity handling
   encoded_matrix = picasso.encode_cnvs_as_ternary(character_matrix)
   print(f'Original: {character_matrix.shape[1]} features')
   print(f'Encoded: {encoded_matrix.shape[1]} features')

**Feature Filtering for Noise Reduction & Performance Improvements**

.. code-block:: python

   # Use encoded data for noisy data handling
   data = encoded_matrix
   
   # Remove features with very low variance (uninformative)
   print(f'Features before filtering: {data.shape[1]}')
   
   # Calculate modal proportion for each feature
   modal_proportions = (data.values == data.mode(axis=0).values).mean(axis=0)
   
   # Keep features where <99% of cells have the modal value
   informative_features = modal_proportions < 0.99
   filtered_data = data.loc[:, informative_features]
   
   print(f'Features after filtering: {filtered_data.shape[1]}')
   print(f'Removed {data.shape[1] - filtered_data.shape[1]} uninformative features')

**Conservative Parameter Settings**

.. code-block:: python

   # Use conservative parameters for noisy data
   model = picasso.Picasso(
       filtered_data,
       min_depth=2,                    # Force minimum depth to explore structure
       max_depth=12,                   # Limit depth to prevent overfitting
       min_clone_size=50,              # Larger clones for noise robustness
       terminate_by='BIC',             # Use conservative BIC-based termination
       bic_penalty_strength=1.2        # Stronger penalty against complexity
   )
   
   print("Fitting model with conservative parameters...")
   model.fit()
   
   # Analyze results
   phylogeny = model.get_phylogeny()
   clone_assignments = model.get_clone_assignments()
   
   print(f"Conservative approach: {len(phylogeny.get_leaves())} clones")
   print("Clone size distribution:")
   print(clone_assignments['clone_id'].value_counts().describe())

Example 3: Advanced Tree Analysis with CloneTree Class
-------------------------------------------------------

This example shows how to extract detailed phylogenetic information from PICASSO results.

**Comprehensive Tree Analysis**

.. code-block:: python

   # Start with a fitted model (from previous examples)
   data = picasso.load_data()
   model = picasso.Picasso(data, min_clone_size=10)
   model.fit()
   
   # Create CloneTree with modal aggregation
   tree = picasso.CloneTree(
       model.get_phylogeny(),
       model.get_clone_assignments(),
       data,
       clone_aggregation='mode'  # Use modal values for clone profiles
   )
   
   # Root the tree
   outgroup = tree.get_most_ancestral_clone()
   tree.root_tree(outgroup)
   print(f"Tree rooted at: {outgroup}")

**Clone Profile Analysis**

We can examine the overall CNA profile that characterizes each clone:

.. code-block:: python

   import numpy as np
   
   # Get modal CNA profiles for each clone
   modal_profiles, modal_frequencies = tree.get_modal_clone_profiles()
   print(f"Modal profiles shape: {modal_profiles.shape}")
   
   # Visualize clone profiles
   plt.figure(figsize=(12, 8))
   sns.clustermap(modal_profiles,
                  cmap='coolwarm', center=0,
                  figsize=(12, 8),
                  cbar_kws={'label': 'Modal Copy Number'},
                  col_cluster=False)
   plt.title('Clone CNA Profiles (Modal Values)')
   plt.show()
   
   # Visualize the frequencies of the modal values to get a sense of how noisy the leaves are
   plt.figure(figsize=(12, 8))
   sns.clustermap(modal_frequencies,
               cmap='Blues', vmin=0,
               cbar_kws={'label': 'Modal Frequency'},
               col_cluster=False)
   plt.title('Clone CNA Profiles (Modal Frequencies)')
   plt.show()

**Sample-Level Phylogeny and Tree Statistics**

.. code-block:: python

   # Create sample-level phylogeny (may be large as it contains all cells)
   print("Creating sample phylogeny...")
   sample_tree = tree.get_sample_phylogeny()
   print(f"Sample tree has {len(sample_tree.get_leaves())} leaves")
   
   # For visualization, we'll work with clone tree
   clone_phylogeny = tree.get_clone_phylogeny()
   
   # Tree statistics
   print(f"Clone tree depth: {clone_phylogeny.get_farthest_leaf()[1]}")
   print(f"Number of internal nodes: {len(clone_phylogeny.get_descendants()) - len(clone_phylogeny.get_leaves())}")

Example 4: iTOL Export for Publication Figures
-----------------------------------------------

iTOL is a sophisticated visualization tool for phylogenies. This example shows how to create publication-ready visualizations using iTOL and helper functions.

**Prepare Phylogeny for iTOL**

.. code-block:: python

   # Prepare data
   data = picasso.load_data()
   model = picasso.Picasso(data, min_clone_size=15)
   model.fit()
   
   tree = picasso.CloneTree(model.get_phylogeny(),
                           model.get_clone_assignments(),
                           data)
   outgroup = tree.get_most_ancestral_clone()
   tree.root_tree(outgroup)
   
   # Get cell-level tree for iTOL (use clone tree if too large)
   cell_tree = tree.get_sample_phylogeny()
   newick_string = cell_tree.write()
   
   # Save tree file for iTOL
   with open('cell_phylogeny.nwk', 'w') as f:
       f.write(newick_string)
   
   print("Saved phylogeny to cell_phylogeny.nwk")

**CNA Heatmap Annotation**

.. code-block:: python

   # Create heatmap annotation showing CNA profiles
   heatmap_annotation = picasso.itol_utils.dataframe_to_itol_heatmap(
       data,
       dataset_label="Copy Number Alterations",
       color_min='#053061',  # Dark blue for deletions
       color_max='#67001f'   # Dark red for amplifications
   )
   
   # Save annotation file
   with open('cna_heatmap.txt', 'w') as f:
       f.write(heatmap_annotation)
   
   print("Saved CNA heatmap annotation to cna_heatmap.txt")
   print("First few lines:")
   print('\\n'.join(heatmap_annotation.split('\\n')[:10]))

**Metadata Color Strips**

.. code-block:: python

   # Create sample metadata for demonstration
   clone_assignments = model.get_clone_assignments()
   
   # Simulate tissue sites
   np.random.seed(42)  # For reproducibility
   sites = np.random.choice(['Primary', 'Metastasis_1', 'Metastasis_2', 'Normal'],
                           size=len(clone_assignments))
   sites_series = pd.Series(sites, index=clone_assignments.index, name='Tissue_Site')
   
   # Define color mapping
   site_colors = {
       'Primary': '#e41a1c',
       'Metastasis_1': '#377eb8',
       'Metastasis_2': '#4daf4a',
       'Normal': '#984ea3'
   }
   
   # Create color strip annotation
   colorstrip_annotation = picasso.itol_utils.dataframe_to_itol_colorstrip(
       sites_series,
       site_colors,
       dataset_label='Tissue Site'
   )
   
   with open('tissue_sites.txt', 'w') as f:
       f.write(colorstrip_annotation)
   
   print("Saved tissue site annotation to tissue_sites.txt")

**Clone Composition Stacked Bars**

.. code-block:: python

   # Analyze tissue composition within each clone
   clone_tissue_data = clone_assignments.merge(sites_series,
                                              left_index=True,
                                              right_index=True)
   
   # Calculate proportions of each tissue type within each clone
   site_proportions = (clone_tissue_data.groupby('clone_id')['Tissue_Site']
                       .value_counts(normalize=True)
                       .unstack(fill_value=0))
   
   print("Tissue proportions by clone:")
   print(site_proportions.head())
   
   # Create stacked bar annotation for clone tree
   stackedbar_annotation = picasso.itol_utils.dataframe_to_itol_stackedbar(
       site_proportions,
       site_colors,
       dataset_label='Tissue Composition'
   )
   
   with open('clone_composition.txt', 'w') as f:
       f.write(stackedbar_annotation)
   
   print("Saved clone composition annotation to clone_composition.txt")

**iTOL Workflow Summary**

.. code-block:: text

   Files created for iTOL:
   1. cell_phylogeny.nwk - Main phylogenetic tree
   2. cna_heatmap.txt - CNA profile heatmap
   3. tissue_sites.txt - Tissue site color strips
   4. clone_composition.txt - Clone composition stacked bars
   
   Steps for iTOL visualization:
   1. Go to https://itol.embl.de/
   2. Upload cell_phylogeny.nwk
   3. Drag and drop annotation files to add visualizations
   4. Customize colors, labels, and layout
   5. Export high-resolution figures
   
   Pro tip: Use clone tree instead of cell tree for large datasets
            to improve iTOL performance and readability

Summary and Best Practices
---------------------------

**Key Takeaways**

1. **Start Simple**: Begin with default parameters and basic workflow
2. **Understand Your Data**: Examine noise levels, feature characteristics, and data quality
3. **Data Preprocessing**: Use ternary encoding and feature filtering for noisy datasets
4. **Parameter Tuning**: Choose conservative parameters for noisy data
5. **Tree Analysis**: Use CloneTree for detailed phylogenetic analysis
6. **Visualization**: Use iTOL for publication-ready phylogenetic figures

**Parameter Selection Guidelines**

Here are general guidelines for parameter selection based on data quality:

================== =============== =============== =================
Parameter          Clean Data      Noisy Data      Very Noisy Data
================== =============== =============== =================
min_clone_size     5-15            15-50           50-100
confidence_thresh  0.7-0.8         0.8-0.85        0.85-0.95
max_depth          unlimited       10-15           8-12
terminate_by       BIC/probability probability     BIC
bic_penalty        1.0             1.0-1.2         1.2-1.5
================== =============== =============== =================

**Common Pitfalls to Avoid**

- Using too small ``min_clone_size`` with noisy data (leads to over-fitting)
- Setting ``max_depth`` too high with noisy data (computational burden, over-fitting)
- Ignoring feature filtering for high-dimensional noisy datasets
- Not using ternary encoding for complex copy number data
- Skipping data quality assessment before parameter selection

**Next Steps**

For more advanced usage, consult the :doc:`api` documentation for detailed parameter descriptions and method specifications.