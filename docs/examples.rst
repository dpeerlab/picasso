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
   sns.heatmap(character_matrix.iloc[:100, :20], 
               cmap='RdBu_r', center=0, 
               cbar_kws={'label': 'Copy Number'})
   plt.title('Copy Number Alterations (first 100 cells, 20 features)')
   plt.xlabel('Genomic Features')
   plt.ylabel('Cells')
   plt.tight_layout()
   plt.show()

**Basic Phylogeny Reconstruction**

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
   print(f"Clone size distribution:")
   print(clone_assignments['clone_id'].value_counts().head())

**Tree Analysis and Visualization**

.. code-block:: python

   # Create CloneTree for advanced analysis
   tree = picasso.CloneTree(phylogeny, clone_assignments, character_matrix)
   
   # Root the tree at the most ancestral clone
   outgroup = tree.get_most_ancestral_clone()
   tree.root_tree(outgroup)
   print(f"Tree rooted at clone: {outgroup}")
   
   # Generate visualizations
   tree.plot_clone_sizes(figsize=(10, 6))
   tree.plot_alterations(figsize=(12, 8))
   
   # Get clone phylogeny as Newick string for external tools
   clone_tree = tree.get_clone_phylogeny()
   print("Newick format (first 100 characters):")
   print(clone_tree.write()[:100] + "...")

Example 2: Handling Very Noisy scRNA-seq Data
----------------------------------------------

This example shows how to handle very noisy CNA data typically obtained from scRNA-seq inference.

**Optional: Ternary Encoding for Complex CNAs**

.. code-block:: python

   # Load data
   character_matrix = picasso.load_data()
   
   # Encode complex CNAs as ternary values for better similarity handling
   encoded_matrix = picasso.encode_cnvs_as_ternary(character_matrix)
   print(f'Original: {character_matrix.shape[1]} features')
   print(f'Encoded: {encoded_matrix.shape[1]} features')
   
   # Visualize encoding difference
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
   
   sns.heatmap(character_matrix.iloc[:50, :10], ax=ax1, 
               cmap='RdBu_r', center=0)
   ax1.set_title('Original Data')
   
   sns.heatmap(encoded_matrix.iloc[:50, :20], ax=ax2, 
               cmap='RdBu_r', center=0)
   ax2.set_title('Ternary Encoded Data')
   
   plt.tight_layout()
   plt.show()

**Feature Filtering for Noise Reduction**

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
       min_clone_size=20,              # Larger clones for noise robustness
       terminate_by='probability',     # Use confidence-based termination
       assignment_confidence_threshold=0.85,  # Higher confidence requirement
       assignment_confidence_proportion=0.95, # Most cells must be confident
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

Example 3: Parameter Tuning and Model Selection
------------------------------------------------

This example demonstrates how to systematically explore parameter space to find optimal settings.

**Systematic Parameter Exploration**

.. code-block:: python

   import numpy as np
   from itertools import product
   
   # Load and prepare data
   data = picasso.load_data()
   
   # Parameter grid
   min_clone_sizes = [5, 10, 20, 50]
   confidence_thresholds = [0.7, 0.8, 0.9]
   
   results = []
   
   print("Exploring parameter combinations...")
   for min_size, threshold in product(min_clone_sizes, confidence_thresholds):
       print(f"Testing min_clone_size={min_size}, threshold={threshold}")
       
       model = picasso.Picasso(
           data,
           min_clone_size=min_size,
           assignment_confidence_threshold=threshold,
           max_depth=10  # Limit for faster exploration
       )
       
       model.fit()
       n_clones = len(model.get_phylogeny().get_leaves())
       
       results.append({
           'min_clone_size': min_size,
           'confidence_threshold': threshold, 
           'n_clones': n_clones
       })
       
       print(f"  -> {n_clones} clones")

**Results Analysis and Visualization**

.. code-block:: python

   # Convert to DataFrame for analysis
   results_df = pd.DataFrame(results)
   
   # Create pivot table for heatmap
   pivot_table = results_df.pivot(
       index='min_clone_size', 
       columns='confidence_threshold', 
       values='n_clones'
   )
   
   # Visualize results
   plt.figure(figsize=(8, 6))
   sns.heatmap(pivot_table, annot=True, fmt='d', 
               cmap='viridis', cbar_kws={'label': 'Number of Clones'})
   plt.title('Parameter Exploration Results')
   plt.ylabel('Minimum Clone Size')
   plt.xlabel('Assignment Confidence Threshold')
   plt.tight_layout()
   plt.show()
   
   # Find optimal parameters (example: moderate number of clones)
   target_clones = 50  # Adjust based on your expectations
   results_df['distance_to_target'] = abs(results_df['n_clones'] - target_clones)
   optimal = results_df.loc[results_df['distance_to_target'].idxmin()]
   
   print("Optimal parameters for ~{} clones:".format(target_clones))
   print(f"  min_clone_size: {optimal['min_clone_size']}")
   print(f"  confidence_threshold: {optimal['confidence_threshold']}")
   print(f"  actual_clones: {optimal['n_clones']}")

Example 4: Advanced Tree Analysis with CloneTree
-------------------------------------------------

This example shows advanced phylogenetic analysis capabilities.

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

.. code-block:: python

   # Get modal CNA profiles for each clone
   modal_profiles = tree.get_modal_clone_profiles()
   print(f"Modal profiles shape: {modal_profiles.shape}")
   
   # Visualize clone profiles
   plt.figure(figsize=(12, 8))
   sns.clustermap(modal_profiles, 
                  cmap='RdBu_r', center=0,
                  figsize=(12, 8),
                  cbar_kws={'label': 'Modal Copy Number'})
   plt.title('Clone CNA Profiles (Modal Values)')
   plt.show()
   
   # Calculate profile similarities
   from scipy.spatial.distance import pdist, squareform
   distances = pdist(modal_profiles, metric='hamming')
   distance_matrix = squareform(distances)
   distance_df = pd.DataFrame(
       distance_matrix, 
       index=modal_profiles.index,
       columns=modal_profiles.index
   )
   
   print("Most similar clones:")
   # Find most similar clone pairs
   mask = np.triu(np.ones_like(distance_matrix, dtype=bool), k=1)
   distance_matrix_masked = distance_matrix.copy()
   distance_matrix_masked[~mask] = np.nan
   
   min_idx = np.nanargmin(distance_matrix_masked)
   i, j = np.unravel_index(min_idx, distance_matrix_masked.shape)
   clone1, clone2 = modal_profiles.index[i], modal_profiles.index[j]
   similarity = 1 - distance_matrix[i, j]
   
   print(f"  {clone1} and {clone2}: {similarity:.3f} similarity")

**Evolutionary Change Analysis**

.. code-block:: python

   # Infer evolutionary changes along tree branches
   changes = tree.infer_evolutionary_changes()
   
   print(f"Detected {len(changes)} evolutionary events")
   print("Sample evolutionary changes:")
   for i, (node, change_dict) in enumerate(changes.items()):
       if i >= 3:  # Show first 3
           break
       print(f"  Node {node}:")
       for feature, change in list(change_dict.items())[:3]:
           print(f"    {feature}: {change}")

**Sample-Level Phylogeny**

.. code-block:: python

   # Create sample-level phylogeny (may be large)
   print("Creating sample phylogeny...")
   sample_tree = tree.get_sample_phylogeny()
   print(f"Sample tree has {len(sample_tree.get_leaves())} leaves")
   
   # For visualization, we'll work with clone tree
   clone_phylogeny = tree.get_clone_phylogeny()
   
   # Tree statistics
   print(f"Clone tree depth: {clone_phylogeny.get_farthest_leaf()[1]}")
   print(f"Number of internal nodes: {len(clone_phylogeny.get_descendants()) - len(clone_phylogeny.get_leaves())}")

Example 5: iTOL Export for Publication Figures
-----------------------------------------------

This example shows how to create publication-ready visualizations using iTOL.

**Basic iTOL Annotations**

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

**iTOL Visualization Workflow**

.. code-block:: python

   print("\\n" + "="*50)
   print("ITOL VISUALIZATION WORKFLOW")
   print("="*50)
   print()
   print("Files created for iTOL:")
   print("1. cell_phylogeny.nwk - Main phylogenetic tree")
   print("2. cna_heatmap.txt - CNA profile heatmap")  
   print("3. tissue_sites.txt - Tissue site color strips")
   print("4. clone_composition.txt - Clone composition stacked bars")
   print()
   print("Steps for iTOL visualization:")
   print("1. Go to https://itol.embl.de/")
   print("2. Upload cell_phylogeny.nwk")
   print("3. Drag and drop annotation files to add visualizations")
   print("4. Customize colors, labels, and layout") 
   print("5. Export high-resolution figures")
   print()
   print("Pro tip: Use clone tree instead of cell tree for large datasets")
   print("         to improve iTOL performance and readability")

Example 6: Custom Preprocessing Pipeline
-----------------------------------------

This example shows advanced preprocessing for challenging datasets.

**Advanced Feature Selection**

.. code-block:: python

   # Load data
   raw_data = picasso.load_data()
   print(f"Raw data: {raw_data.shape}")
   
   # Step 1: Remove features with extreme modal proportions
   modal_props = (raw_data.values == raw_data.mode(axis=0).values).mean(axis=0)
   
   # Keep features with 5-95% modality (adjustable)
   informative_mask = (modal_props >= 0.05) & (modal_props <= 0.95)
   data_step1 = raw_data.loc[:, informative_mask]
   print(f"After modality filtering: {data_step1.shape}")
   
   # Step 2: Remove features with very low variance
   feature_variance = data_step1.var()
   variance_threshold = feature_variance.quantile(0.1)  # Bottom 10%
   
   high_var_mask = feature_variance > variance_threshold
   data_step2 = data_step1.loc[:, high_var_mask]
   print(f"After variance filtering: {data_step2.shape}")
   
   # Step 3: Optional - correlation-based feature reduction
   corr_matrix = data_step2.corr().abs()
   
   # Find highly correlated feature pairs
   high_corr_pairs = []
   for i in range(len(corr_matrix.columns)):
       for j in range(i+1, len(corr_matrix.columns)):
           if corr_matrix.iloc[i, j] > 0.95:  # Highly correlated
               high_corr_pairs.append((corr_matrix.columns[i], corr_matrix.columns[j]))
   
   print(f"Found {len(high_corr_pairs)} highly correlated feature pairs")
   
   # Remove one feature from each highly correlated pair
   features_to_remove = set()
   for feat1, feat2 in high_corr_pairs:
       # Keep the feature with higher variance
       if data_step2[feat1].var() < data_step2[feat2].var():
           features_to_remove.add(feat1)
       else:
           features_to_remove.add(feat2)
   
   final_features = [f for f in data_step2.columns if f not in features_to_remove]
   processed_data = data_step2[final_features]
   print(f"Final processed data: {processed_data.shape}")

**Quality Control Analysis**

.. code-block:: python

   # Analyze preprocessing impact
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   
   # Feature count by processing step
   steps = ['Raw', 'Modality Filter', 'Variance Filter', 'Correlation Filter']
   counts = [raw_data.shape[1], data_step1.shape[1], 
            data_step2.shape[1], processed_data.shape[1]]
   
   axes[0,0].bar(steps, counts, color=['red', 'orange', 'yellow', 'green'])
   axes[0,0].set_title('Features Retained by Processing Step')
   axes[0,0].set_ylabel('Number of Features')
   plt.setp(axes[0,0].xaxis.get_majorticklabels(), rotation=45)
   
   # Variance distributions
   axes[0,1].hist(raw_data.var(), bins=50, alpha=0.7, label='Raw', color='red')
   axes[0,1].hist(processed_data.var(), bins=50, alpha=0.7, label='Processed', color='green')
   axes[0,1].set_title('Feature Variance Distribution')
   axes[0,1].set_xlabel('Variance')
   axes[0,1].legend()
   
   # Modal proportion distributions  
   raw_modal = (raw_data.values == raw_data.mode(axis=0).values).mean(axis=0)
   proc_modal = (processed_data.values == processed_data.mode(axis=0).values).mean(axis=0)
   
   axes[1,0].hist(raw_modal, bins=50, alpha=0.7, label='Raw', color='red')
   axes[1,0].hist(proc_modal, bins=50, alpha=0.7, label='Processed', color='green')
   axes[1,0].set_title('Modal Proportion Distribution')
   axes[1,0].set_xlabel('Modal Proportion')
   axes[1,0].legend()
   
   # Data range comparison
   axes[1,1].boxplot([raw_data.values.flatten(), processed_data.values.flatten()],
                    labels=['Raw', 'Processed'])
   axes[1,1].set_title('Copy Number Value Distribution')
   axes[1,1].set_ylabel('Copy Number')
   
   plt.tight_layout()
   plt.show()

**Optimized Model Fitting**

.. code-block:: python

   # Fit model with processed data
   print("Fitting PICASSO with processed data...")
   
   model_processed = picasso.Picasso(
       processed_data,
       min_clone_size=15,
       assignment_confidence_threshold=0.8,
       assignment_confidence_proportion=0.85
   )
   
   model_processed.fit()
   
   # Compare with raw data model
   print("Fitting PICASSO with raw data for comparison...")
   model_raw = picasso.Picasso(
       raw_data,
       min_clone_size=15,
       assignment_confidence_threshold=0.8,
       assignment_confidence_proportion=0.85,
       max_depth=model_processed.depth  # Match depth for fair comparison
   )
   
   model_raw.fit()
   
   # Compare results
   print("\\nComparison Results:")
   print(f"Raw data: {len(model_raw.get_phylogeny().get_leaves())} clones")
   print(f"Processed data: {len(model_processed.get_phylogeny().get_leaves())} clones")
   
   # Analyze clone size distributions
   raw_assignments = model_raw.get_clone_assignments()
   proc_assignments = model_processed.get_clone_assignments()
   
   raw_sizes = raw_assignments['clone_id'].value_counts()
   proc_sizes = proc_assignments['clone_id'].value_counts()
   
   print(f"\\nClone size statistics:")
   print(f"Raw - Mean: {raw_sizes.mean():.1f}, Std: {raw_sizes.std():.1f}")
   print(f"Processed - Mean: {proc_sizes.mean():.1f}, Std: {proc_sizes.std():.1f}")

Summary and Best Practices
---------------------------

**Key Takeaways**

1. **Start Simple**: Begin with default parameters and basic workflow
2. **Understand Your Data**: Examine noise levels, feature characteristics, and data quality
3. **Parameter Tuning**: Systematically explore parameter space for optimal results
4. **Preprocessing**: Filter uninformative features for noisy datasets
5. **Validation**: Compare results across different parameter settings
6. **Visualization**: Use iTOL for publication-ready phylogenetic figures

**Parameter Selection Guidelines**

Here are some rough guidelines for parameter selection. These are not strict rules and should be used as a starting point;
the best parameters will depend on the specific dataset and the user's goals.

================== =============== =============== =================
Parameter          Clean Data      Noisy Data      Very Noisy Data
================== =============== =============== =================
min_clone_size     10-50           50-100          100-200
confidence_thresh  0.7-0.8         0.8-0.85        0.85-0.95
max_depth          unlimited       10-15           8-12
terminate_by       BIC/probability probability     probability
bic_penalty        1.0             1.0-1.2         1.2-1.5
================== =============== =============== =================

**Pitfalls to Avoid**

- Using too small ``min_clone_size`` with noisy data (leads to over-fitting)
- Setting ``max_depth`` too high with noisy data (computational burden, over-fitting)
- Ignoring feature filtering for high-dimensional noisy datasets
- Not validating results across different parameter settings

**Next Steps**

For more advanced usage, consult the :doc:`api` documentation for detailed parameter descriptions and method specifications.