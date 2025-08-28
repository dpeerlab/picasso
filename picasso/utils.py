"""
Utilities: Data preprocessing and loading functions for PICASSO.

This module provides utility functions for preprocessing copy number alteration (CNA)
data and loading example datasets. It includes specialized functions for handling
noisy scRNA-seq-inferred CNA data and converting complex copy number states into
formats suitable for phylogenetic analysis.

Functions
---------
encode_cnvs_as_ternary
    Convert integer CNA data to ternary encoding for improved phylogenetic inference.
load_data  
    Load example CNA dataset for testing and demonstration purposes.

Examples
--------
Data preprocessing workflow:

>>> from picasso import Picasso, load_data, encode_cnvs_as_ternary
>>> 
>>> # Load example dataset
>>> cna_data = load_data()
>>> print(f"Loaded data: {cna_data.shape}")
>>> 
>>> # Optional: Convert to ternary encoding for complex copy number states
>>> ternary_data = encode_cnvs_as_ternary(cna_data)
>>> print(f"Ternary encoded: {ternary_data.shape}")
>>> 
>>> # Use with PICASSO
>>> picasso = Picasso(cna_data, min_clone_size=8)
>>> picasso.fit()

Notes
-----
These utilities are specifically designed for:
- Handling noisy scRNA-seq-inferred CNA data
- Converting complex copy number states to phylogeny-compatible formats
- Providing realistic example data for algorithm development
- Supporting data preprocessing workflows

See Also
--------
Picasso : Main phylogenetic inference class
CloneTree : Analysis and visualization of phylogenetic results
"""

import numpy as np
import pandas as pd
import os
from typing import Union

def encode_cnvs_as_ternary(data: Union[pd.DataFrame, np.ndarray]) -> pd.DataFrame:
    """
    Convert CNA data to ternary encoding for phylogenetic analysis.
    
    Transforms integer copy number alteration (CNA) data into a ternary format
    suitable for phylogenetic inference algorithms like PICASSO. This encoding 
    is particularly useful for handling complex copy number states and ensuring
    compatibility with categorical mixture models.

    Parameters
    ----------
    data : pd.DataFrame or np.ndarray
        Input CNA data where rows represent cells/samples and columns represent
        genomic features. Values should be integers representing copy number states
        (e.g., 0=deletion, 1=neutral, 2=single amplification, 3=double amplification).
        Can handle both positive and negative copy number values. 

    Returns
    -------
    pd.DataFrame
        Ternary-encoded DataFrame with values in {-1, 0, 1}. The number of columns
        is expanded based on the maximum absolute value in each original column.
        Column names follow the pattern 'original_column-position' (e.g., 'chr1p-1', 'chr1p-2').

    Examples
    --------
    Basic encoding of copy number states:
    
    >>> import pandas as pd
    >>> import numpy as np
    >>> from picasso.utils import encode_cnvs_as_ternary
    >>> 
    >>> # Create sample CNA data
    >>> cna_data = pd.DataFrame({
    ...     'chr1p': [0, 1, 2, 3],
    ...     'chr2q': [0, 0, 1, 2]
    ... }, index=['Cell_A', 'Cell_B', 'Cell_C', 'Cell_D'])
    >>> 
    >>> print(cna_data)
           chr1p  chr2q
    Cell_A     0     0
    Cell_B     -1     0
    Cell_C     2     1
    Cell_D     3     2
    
    >>> # Encode to ternary format
    >>> ternary_data = encode_cnvs_as_ternary(cna_data)
    >>> print(ternary_data)
           chr1p-1  chr1p-2  chr1p-3  chr2q-1  chr2q-2
    Cell_A       0        0        0        0        0
    Cell_B      -1        0        0        0        0
    Cell_C       1        1        0        1        0
    Cell_D       1        1        1        1        1

    Handling deletions (negative values):
    
    >>> # Data with deletions
    >>> cna_with_dels = pd.DataFrame({
    ...     'chr3p': [-2, -1, 0, 1, 2],
    ... }, index=[f'Cell_{i}' for i in range(5)])
    >>> 
    >>> ternary_dels = encode_cnvs_as_ternary(cna_with_dels)
    >>> print(ternary_dels)
           chr3p-1  chr3p-2
    Cell_0      -1       -1
    Cell_1      -1        0
    Cell_2       0        0
    Cell_3       1        0
    Cell_4       1        1

    Notes
    -----
    **Encoding Rules**:
    - Positive integers n are encoded as n ones followed by zeros: [1, 1, ..., 1, 0, 0, ...]
    - Negative integers -n are encoded as n negative ones: [-1, -1, ..., -1]
    - Zero values are encoded as all zeros: [0, 0, ...]
    - Column width is determined by the maximum absolute value in each original column
    
    **Use Cases**:
    - Preprocessing CNA data for PICASSO phylogenetic inference
    - Converting complex copy number states to categorical format
    - Ensuring proper handling of amplifications and deletions in mixture models
    
    **Performance Considerations**:
    - Output size scales with maximum copy number values
    - Memory usage increases significantly for high-amplitude CNAs
    - Consider binning extreme values before encoding for very noisy data. We recommend binning into 'amplified' and 'highly amplified' categories.

    Raises
    ------
    ValueError
        If input data cannot be converted to integer format.

    See Also
    --------
    Picasso : Main phylogenetic inference class that accepts ternary-encoded data
    load_data : Function to load example CNA datasets
    """

    # If input is a numpy array, convert it to a DataFrame
    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data)

    data = data.astype(int)

    # Initialize a list to hold the binary encoded columns
    binary_encoded_cols = []
    column_names = []

    # Process each column independently
    for col in data.columns:
        col_data = data[col]
        # Get the maximum magnitude in the column
        max_val = np.max(np.abs(col_data))

        # Initialize an empty list to hold the binary encoded values for the column
        binary_col = []

        # Encode each value in the column
        for val in col_data:
            if val >= 0:
                binary_val = [1] * val + [0] * (max_val - val)
            else:
                binary_val = [-1] * abs(val)
            binary_col.append(binary_val)

        # Determine the length needed for padding
        max_length = max(len(b) for b in binary_col)

        # Pad binary_col to ensure uniform length
        padded_col = [np.pad(b, (0, max_length - len(b)), 'constant') for b in binary_col]

        # Convert the padded column to a numpy array
        padded_col = np.array(padded_col)

        # Add the binary encoded columns to the list
        for i in range(max_length):
            binary_encoded_cols.append(padded_col[:, i])
            column_names.append(f"{col}-{i + 1}")

    # Combine all binary encoded columns into a DataFrame
    binary_encoded_df = pd.DataFrame(np.column_stack(binary_encoded_cols), columns=column_names)
    binary_encoded_df.index = data.index

    return binary_encoded_df


def load_data() -> pd.DataFrame:
    """
    Load example single-cell copy number alteration (CNA) dataset.
    
    Provides a sample dataset of inferred CNAs from single-cell RNA sequencing data
    for testing and demonstration purposes. This dataset represents the type of noisy,
    inferred CNA data that PICASSO is designed to handle.

    Returns
    -------
    pd.DataFrame
        Example CNA dataset with cells as rows and genomic features as columns.
        Values represent inferred copy number states, typically integers where:
        - 0 indicates deletions/loss
        - 1 indicates neutral copy number  
        - 2+ indicates amplifications/gains
        Index contains cell/sample identifiers, columns contain feature names.

    Examples
    --------
    Load and explore the example dataset:
    
    >>> from picasso import Picasso, load_data
    >>> 
    >>> # Load example data
    >>> cna_data = load_data()
    >>> print(f"Dataset shape: {cna_data.shape}")
    >>> print(f"Copy number range: {cna_data.min().min()} to {cna_data.max().max()}")
    >>> print("First few rows:")
    >>> print(cna_data.head())
    >>> 
    >>> # Use with PICASSO
    >>> picasso = Picasso(cna_data, min_clone_size=5)
    >>> picasso.fit()

    Inspect data characteristics:
    
    >>> # Check for missing values
    >>> print(f"Missing values: {cna_data.isnull().sum().sum()}")
    >>> 
    >>> # Distribution of copy number states
    >>> print("Copy number state distribution:")
    >>> print(cna_data.values.flatten().astype(int))
    >>> 
    >>> # Feature-wise statistics
    >>> print("Per-feature statistics:")
    >>> print(cna_data.describe())

    Notes
    -----
    **Dataset Characteristics**:
    - Representative of scRNA-seq-inferred CNA data
    - Contains typical noise patterns and artifacts  
    - Suitable for algorithm testing and parameter tuning
    - May include both amplifications and deletions
    
    **Data Origin**:
    - Loaded from sample_data/cnas.txt in the package directory
    - Tab-separated format with sample IDs as first column
    - Preprocessed to remove extreme outliers and artifacts
    
    **Intended Use**:
    - Algorithm development and testing
    - Parameter optimization for noisy datasets
    - Tutorial and documentation examples
    - Benchmarking against other methods

    Raises
    ------
    FileNotFoundError
        If the sample data file cannot be located in the expected directory.
    pd.errors.EmptyDataError
        If the data file is empty or corrupted.

    See Also
    --------
    Picasso : Main phylogenetic inference class for analyzing the loaded data
    encode_cnvs_as_ternary : Preprocessing function for complex copy number states
    CloneTree : Class for visualizing and analyzing phylogenetic results
    """

    # Load the example dataset
    # Get path to sample data within the package
    package_dir = os.path.dirname(os.path.abspath(__file__))
    data = pd.read_csv(f'{package_dir}/sample_data/cnas.txt', sep='\t', index_col=0)
    return data


# Define public API
__all__ = ['encode_cnvs_as_ternary', 'load_data']
