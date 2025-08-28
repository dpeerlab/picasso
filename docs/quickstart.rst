Quick Start Guide
=================

This guide will get you up and running with PICASSO in just a few minutes. PICASSO reconstructs phylogenetic trees from noisy copy number alteration (CNA) data derived from single-cell RNA sequencing.

.. note::
   Before starting, make sure you have installed PICASSO following the :doc:`installation` guide.

5-Minute Tutorial
-----------------

**Step 1: Load Your Data**

PICASSO works with CNA data as a pandas DataFrame where rows are cells and columns are genomic features:

.. code-block:: python

   import picasso
   import pandas as pd
   
   # Load example data (or use your own)
   cna_data = picasso.load_data()
   print(f"Dataset: {cna_data.shape[0]} cells × {cna_data.shape[1]} features")
   
   # Examine the data structure
   print(cna_data.head())

The data should contain integer copy number states (e.g., -1=deletion, 0=neutral, 1=amplification).

**Step 2: Basic Phylogeny Reconstruction**

For most datasets, the default parameters work well:

.. code-block:: python

   # Initialize PICASSO with default parameters
   model = picasso.Picasso(cna_data)
   
   # Reconstruct the phylogeny
   model.fit()
   
   # Extract results
   phylogeny = model.get_phylogeny()
   clone_assignments = model.get_clone_assignments()
   
   print(f"Reconstructed {len(phylogeny.get_leaves())} clones")

**Step 3: Analyze Results**

Use CloneTree for enhanced analysis and visualization:

.. code-block:: python

   # Create tree analyzer
   tree = picasso.CloneTree(phylogeny, clone_assignments, cna_data)
   
   # Root the tree at the most ancestral clone
   outgroup = tree.get_most_ancestral_clone()
   tree.root_tree(outgroup)
   
   # Visualize clone sizes and alterations
   tree.plot_clone_sizes()
   tree.plot_alterations()

**Step 4: Export Results**

Export for further analysis or publication:

.. code-block:: python

   # Get clone phylogeny as Newick string
   clone_tree = tree.get_clone_phylogeny()
   newick_string = clone_tree.write()
   
   # Save results
   clone_assignments.to_csv('clone_assignments.csv')
   
   # Export for iTOL visualization
   heatmap_annotation = picasso.itol_utils.dataframe_to_itol_heatmap(cna_data)
   with open('itol_heatmap.txt', 'w') as f:
       f.write(heatmap_annotation)

That's it! You've successfully reconstructed a phylogeny from CNA data.

Understanding Your Data
-----------------------

**Input Format**

PICASSO expects a pandas DataFrame with:
- **Rows**: Individual cells/samples
- **Columns**: Genomic features (chromosome arms, genes, bins)
- **Values**: Integer copy number states

Common encodings:
- ``-2, -1``: Deletions (homozygous, heterozygous)
- ``0``: Neutral copy number
- ``1, 2, 3+``: Amplifications (single, double, triple+)

**Data Quality Considerations**

PICASSO is designed for noisy scRNA-seq-inferred CNAs, but data quality affects results:

.. code-block:: python

   # Check data characteristics
   print(f"Copy number range: {cna_data.min().min()} to {cna_data.max().max()}")
   print(f"Missing values: {cna_data.isnull().sum().sum()}")
   print(f"Feature variance: {cna_data.var().describe()}")

**Feature Filtering (Optional)**

Remove uninformative features to speed up analysis:

.. code-block:: python

   # Remove features with >99% modal values (optional)
   n_features_before = cna_data.shape[1]
   modal_threshold = 0.99
   
   feature_modality = (cna_data.values == cna_data.mode(axis=0).values).mean(axis=0)
   informative_features = feature_modality < modal_threshold
   
   cna_filtered = cna_data.loc[:, informative_features]
   print(f"Filtered from {n_features_before} to {cna_filtered.shape[1]} features")

Parameter Selection
-------------------

**For Clean Data**

If your CNA data has low noise (e.g., from scDNA-seq or well-validated inference):

.. code-block:: python

   model = picasso.Picasso(
       cna_data,
       min_clone_size=5,              # Smaller clones OK
       assignment_confidence_threshold=0.7,  # Lower confidence OK
       terminate_by='BIC'             # BIC-based termination
   )

**For Noisy Data** 

If your data is very noisy (typical for scRNA-seq-inferred CNAs):

.. code-block:: python

   model = picasso.Picasso(
       cna_data,
       min_clone_size=20,             # Larger clones for robustness  
       max_depth=10,                  # Limit depth to avoid over-fitting
       assignment_confidence_threshold=0.85,  # Higher confidence required
       assignment_confidence_proportion=0.9,  # Most cells must be confident
       terminate_by='probability'     # Confidence-based termination
   )

**Parameter Guidelines**

- **min_clone_size**: 5-10 for clean data, 10-50+ for noisy data
- **assignment_confidence_threshold**: 0.7 for clean data, 0.8-0.9 for noisy data
- **max_depth**: Unlimited for clean data, 8-15 for noisy data to prevent over-fitting
- **terminate_by**: 'BIC' for clean data, 'probability' for noisy data

Common Workflows
----------------

**Workflow 1: Standard Analysis**

.. code-block:: python

   # 1. Load and examine data
   data = picasso.load_data()
   
   # 2. Reconstruct phylogeny
   model = picasso.Picasso(data, min_clone_size=10)
   model.fit()
   
   # 3. Analyze results
   tree = picasso.CloneTree(model.get_phylogeny(), 
                           model.get_clone_assignments(), 
                           data)
   
   # 4. Generate visualizations
   tree.plot_clone_sizes()
   tree.plot_alterations()

**Workflow 2: Parameter Exploration**

.. code-block:: python

   # Try different minimum clone sizes
   clone_sizes = [5, 10, 20, 50]
   results = {}
   
   for size in clone_sizes:
       model = picasso.Picasso(data, min_clone_size=size)
       model.fit()
       n_clones = len(model.get_phylogeny().get_leaves())
       results[size] = n_clones
       print(f"Min clone size {size}: {n_clones} clones")

**Workflow 3: Noisy Data Pipeline**

.. code-block:: python

   # 1. Filter features
   modality = (data.values == data.mode(axis=0).values).mean(axis=0)
   filtered_data = data.loc[:, modality < 0.95]
   
   # 2. Use conservative parameters  
   model = picasso.Picasso(
       filtered_data,
       min_clone_size=25,
       max_depth=12,
       assignment_confidence_threshold=0.85,
       assignment_confidence_proportion=0.9
   )
   
   # 3. Fit and analyze
   model.fit()
   tree = picasso.CloneTree(model.get_phylogeny(), 
                           model.get_clone_assignments(),
                           filtered_data)

Output Interpretation
---------------------

**Clone Assignments**

The clone assignments DataFrame shows which clone each cell belongs to:

.. code-block:: python

   assignments = model.get_clone_assignments()
   print(assignments.head())
   
   # Clone size distribution
   clone_sizes = assignments['clone_id'].value_counts()
   print("Clone sizes:", clone_sizes.head())

**Phylogenetic Tree**

The phylogeny represents evolutionary relationships:

.. code-block:: python

   phylogeny = model.get_phylogeny()
   
   # Tree statistics
   print(f"Number of leaves (clones): {len(phylogeny.get_leaves())}")
   print(f"Tree depth: {phylogeny.get_farthest_leaf()[1]}")
   
   # Leaf names correspond to clone IDs
   leaf_names = phylogeny.get_leaf_names()
   print("Clone IDs:", leaf_names[:5])

**Tree Analysis**

CloneTree provides additional insights:

.. code-block:: python

   # Get modal CNA profiles for each clone
   modal_profiles = tree.get_modal_clone_profiles()
   
   # Infer evolutionary changes along branches
   changes = tree.infer_evolutionary_changes()
   print(f"Detected {len(changes)} evolutionary events")

Next Steps
----------

Now that you understand the basics:

1. **Explore Examples**: See :doc:`examples` for detailed tutorials on specific use cases
2. **Parameter Tuning**: Learn how to optimize parameters for your specific data type
3. **Advanced Analysis**: Discover CloneTree's advanced features for phylogenetic analysis
4. **Visualization**: Create publication-ready figures with iTOL integration
5. **API Reference**: Consult :doc:`api` for complete function documentation

**Need Help?**

- Check the :doc:`examples` for similar use cases
- Consult the :doc:`api` for detailed parameter descriptions  
- Visit our `GitHub Issues <https://github.com/yourusername/picasso/issues>`_ for questions
- Review the original publication for algorithmic details