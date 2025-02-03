import numpy as np
import pandas as pd
import os

def encode_cnvs_as_ternary(data):
    """
    Convert an integer pandas DataFrame or numpy array to a ternary DataFrame
    containing only -1, 0, and 1 values. Integer values are encoded as follows:
    - Positive integers are encoded as sequences of 1s followed by 0s
    e.g. If 3 is the maximum (absolute) value in the column, 2 will be encoded as
    [1, 1, 0] and 3 will be encoded as [1, 1, 1]
    - Negative integers are encoded as sequences of -1s followed by 0s
    e.g. If 3 is the maximum (absolute) value in the column, -2 will be encoded as
    [-1, -1, 0] and -3 will be encoded as [-1, -1, -1].

    Parameters:
    data (pd.DataFrame or np.ndarray): Input data.

    Returns:
    pd.DataFrame: Binary encoded DataFrame with appropriate column names.
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


def load_data():
    """
    Load the example dataset containing single-cell copy number variation (CNV) data.

    Returns:
    pd.DataFrame: Example CNV dataset.
    """

    # Load the example dataset
    # Get absolute path to parent directory
    parent_dir = '/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
    data = pd.read_csv(f'{parent_dir}/sample_data/cnas.txt', sep='\t', index_col=0)
    return data
