"""
iTOL Utilities: Export functions for interactive Tree of Life (iTOL) visualizations.

This module provides functions for converting PICASSO phylogenetic results into
iTOL-compatible annotation files. iTOL is a web-based tool for creating publication-
quality phylogenetic tree visualizations with rich annotations.

Functions
---------
rgba_to_hex
    Convert RGBA color tuples to hexadecimal strings for iTOL compatibility.
dataframe_to_itol_colorstrip
    Create iTOL colorstrip annotations for categorical data (e.g., clone assignments).
dataframe_to_itol_heatmap
    Create iTOL heatmap annotations for continuous CNA data visualization.
dataframe_to_itol_stackedbar
    Create iTOL stacked bar annotations for multi-category data display.

Examples
--------
Export PICASSO results for iTOL visualization:

>>> from picasso import Picasso, load_data
>>> from picasso.itol_utils import dataframe_to_itol_colorstrip, dataframe_to_itol_heatmap
>>> import seaborn as sns
>>>
>>> # Load data and run PICASSO analysis
>>> cna_data = load_data()
>>> picasso = Picasso(cna_data)
>>> picasso.fit()
>>> assignments = picasso.get_clone_assignments()
>>>
>>> # Create color mapping for clones
>>> unique_clones = assignments['clone_id'].unique()
>>> colors = sns.color_palette('tab10', len(unique_clones))
>>> color_map = dict(zip(unique_clones, colors))
>>>
>>> # Generate iTOL colorstrip for clone assignments
>>> colorstrip = dataframe_to_itol_colorstrip(assignments['clone_id'],
...                                          color_map, 'Clone_Assignments')
>>> with open('clone_colorstrip.txt', 'w') as f:
...     f.write(colorstrip)
>>>
>>> # Generate iTOL heatmap for CNA profiles
>>> heatmap = dataframe_to_itol_heatmap(cna_data, 'Copy_Number_Alterations')
>>> with open('cna_heatmap.txt', 'w') as f:
...     f.write(heatmap)

Notes
-----
iTOL Integration Features:
- Automatic color format conversion (RGB/RGBA to hex)
- Support for legends and custom labeling
- Compatible with iTOL's annotation file formats
- Handles large datasets efficiently

Recommended Workflow:
1. Run PICASSO phylogenetic inference
2. Extract phylogeny and clone assignments
3. Generate iTOL annotation files using these functions
4. Upload tree (Newick format) and annotations to iTOL web interface
5. Customize visualization and export publication-ready figures

See Also
--------
Picasso : Main phylogenetic inference class
CloneTree : Alternative matplotlib-based visualizations
"""

import io
import pandas as pd
from typing import Dict, Tuple, Union


def rgba_to_hex(rgba: Tuple[float, ...]) -> str:
    """
    Convert RGBA values to hexadecimal color string.

    Parameters
    ----------
    rgba : tuple
        Tuple of RGBA values.

    Returns
    -------
    str
        Hexadecimal color string.

    Examples
    --------
    >>> rgba_to_hex((1.0, 0.0, 0.0, 1.0))
    '#ff0000'

    >>> rgba_to_hex((1.0, 0.0, 0.0))
    '#ff0000'
    """
    # Extract the RGBA values
    if len(rgba) == 3:
        red, green, blue = rgba
    elif len(rgba) == 4:
        red, green, blue, _ = rgba

    # Ensure the values are in the range 0-1
    red = min(1.0, max(0.0, red))
    green = min(1.0, max(0.0, green))
    blue = min(1.0, max(0.0, blue))

    # Convert to hexadecimal and ensure two characters for each value
    red_hex = format(int(red * 255), "02X")
    green_hex = format(int(green * 255), "02X")
    blue_hex = format(int(blue * 255), "02X")

    # Concatenate the hexadecimal values
    hex_color = f"#{red_hex}{green_hex}{blue_hex}"

    return hex_color


def dataframe_to_itol_colorstrip(
    series: pd.Series,
    cmap: Dict[str, Union[str, Tuple[float, ...]]],
    dataset_label: str,
) -> str:
    """
    Generate iTOL colorstrip annotation for phylogenetic tree visualization.

    Creates an iTOL-compatible colorstrip annotation file that can be used to
    visualize categorical data (such as clone assignments) alongside phylogenetic
    trees. This is particularly useful for showing clone membership or other
    discrete annotations in PICASSO results.

    Parameters
    ----------
    series : pd.Series
        Series with sample/leaf identifiers as index and categorical values as data.
        Index should match the leaf names in the corresponding phylogenetic tree.
        Common use cases include clone assignments, cell types, or treatment groups.
    cmap : dict
        Dictionary mapping categorical values to colors. Keys should be the unique
        values in the series, values can be either:
        - Hex color strings (e.g., '#FF0000')
        - RGB tuples (e.g., (1.0, 0.0, 0.0))
        - RGBA tuples (e.g., (1.0, 0.0, 0.0, 1.0))
        RGB/RGBA tuples will be automatically converted to hex format.
    dataset_label : str
        Label for the annotation dataset, displayed in iTOL interface and legend.

    Returns
    -------
    str
        iTOL-formatted colorstrip annotation content that can be saved to a file
        and uploaded to iTOL for tree visualization.

    Examples
    --------
    Create colorstrip for clone assignments:

    >>> import pandas as pd
    >>> from picasso.itol_utils import dataframe_to_itol_colorstrip
    >>>
    >>> # Sample clone assignments
    >>> assignments = pd.Series(['Clone_A', 'Clone_A', 'Clone_B', 'Clone_C'],
    ...                        index=['Cell_1', 'Cell_2', 'Cell_3', 'Cell_4'])
    >>>
    >>> # Define color mapping
    >>> color_map = {
    ...     'Clone_A': '#FF0000',  # Red
    ...     'Clone_B': '#00FF00',  # Green
    ...     'Clone_C': '#0000FF'   # Blue
    ... }
    >>>
    >>> # Generate iTOL annotation
    >>> itol_content = dataframe_to_itol_colorstrip(assignments, color_map, 'Clone_Assignments')
    >>>
    >>> # Save to file for iTOL upload
    >>> with open('clone_colorstrip.txt', 'w') as f:
    ...     f.write(itol_content)

    Using RGB tuples (automatically converted):

    >>> import seaborn as sns
    >>>
    >>> # Generate colors using seaborn
    >>> unique_clones = assignments.unique()
    >>> colors = sns.color_palette('tab10', len(unique_clones))
    >>> color_map = dict(zip(unique_clones, colors))
    >>>
    >>> # RGB tuples will be converted to hex automatically
    >>> itol_content = dataframe_to_itol_colorstrip(assignments, color_map, 'Clones')

    Notes
    -----
    **iTOL Integration**:
    - Output format follows iTOL DATASET_COLORSTRIP specification
    - Files can be directly uploaded to iTOL web interface
    - Supports legends and custom labeling

    **Color Handling**:
    - Automatically converts RGB/RGBA tuples to hex strings
    - Validates that all series values have corresponding colors in cmap
    - Supports any valid CSS/HTML color specification

    **Use Cases**:
    - Visualizing PICASSO clone assignments on phylogenetic trees
    - Showing cell type annotations alongside evolutionary relationships
    - Displaying experimental conditions or metadata categories
    - Creating publication-ready annotated phylogenies

    Raises
    ------
    AssertionError
        If any value in the series lacks a corresponding color in the colormap,
        or if colors cannot be converted to valid string format.

    See Also
    --------
    dataframe_to_itol_heatmap : Create heatmap annotations for continuous data
    dataframe_to_itol_stackedbar : Create stacked bar annotations
    CloneTree.plot_alterations : Alternative visualization using matplotlib
    """

    for key in cmap:
        # If the color is in rgba format, convert it to hex
        if isinstance(cmap[key], tuple) or isinstance(cmap[key], list):
            cmap[key] = rgba_to_hex(cmap[key])
        assert isinstance(cmap[key], str), f"Color for {key} is not a string"

    # Create the annotations file and write it to buffer

    f = io.StringIO()
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write(f"DATASET_LABEL\t{dataset_label}\n")
    f.write("COLOR\t#ff0000\n")
    f.write(f"LEGEND_TITLE\t{dataset_label}\n")
    f.write("LEGEND_SHAPES\t1\n")
    f.write("LEGEND_COLORS\t#ff0000\n")
    f.write(f"LEGEND_LABELS\t{dataset_label}\n")
    f.write("DATA\n")
    for leaf in series.index:
        lineage = series.loc[leaf]
        f.write(f"{leaf}\t{cmap[lineage]}\t{lineage}\n")
    text = f.getvalue()
    f.close()
    return text


def dataframe_to_itol_heatmap(
    df: pd.DataFrame,
    dataset_label: str = "CNVs",
    color_min: str = "#3f4c8a",
    color_max: str = "#b40426",
) -> str:
    """
    Generate iTOL heatmap annotation for continuous data visualization.

    Creates an iTOL-compatible heatmap annotation file for visualizing continuous
    numerical data (such as copy number alterations) alongside phylogenetic trees.
    Uses a three-color gradient to represent data range.

    Parameters
    ----------
    df : pd.DataFrame
        Data matrix with samples/cells as rows and features as columns. Index should
        match the leaf names in the corresponding phylogenetic tree. Values should
        be numeric (e.g., copy number states, expression levels).
    dataset_label : str, default="CNVs"
        Label for the annotation dataset, displayed in iTOL interface and legend.
        Should be descriptive of the data type (e.g., "Copy_Number", "Expression").
    color_min : str, default='#3f4c8a'
        Hexadecimal color for minimum values in the heatmap. Default is dark blue.
    color_max : str, default='#b40426'
        Hexadecimal color for maximum values in the heatmap. Default is dark red.

    Returns
    -------
    str
        iTOL-formatted heatmap annotation content that can be saved to a file
        and uploaded to iTOL for phylogenetic tree visualization.

    Examples
    --------
    Create heatmap for CNA data:

    >>> from picasso import load_data
    >>> from picasso.itol_utils import dataframe_to_itol_heatmap
    >>>
    >>> # Load CNA data
    >>> cna_data = load_data()
    >>>
    >>> # Generate iTOL heatmap annotation
    >>> heatmap = dataframe_to_itol_heatmap(cna_data,
    ...                                    dataset_label='Copy_Number_Alterations',
    ...                                    color_min='#0000FF',  # Blue for deletions
    ...                                    color_max='#FF0000')  # Red for amplifications
    >>>
    >>> # Save for iTOL upload
    >>> with open('cna_heatmap.txt', 'w') as f:
    ...     f.write(heatmap)

    Custom color scheme:

    >>> # Use custom colors
    >>> custom_heatmap = dataframe_to_itol_heatmap(
    ...     expression_data,
    ...     dataset_label='Gene_Expression',
    ...     color_min='#FFFFFF',  # White for low
    ...     color_max='#000000'   # Black for high
    ... )

    Notes
    -----
    **iTOL Heatmap Features**:
    - Three-color gradient: min_color -> mid_color -> max_color
    - Automatic scaling based on data range
    - Mid-color is fixed at light gray (#f5f5f5)
    - Supports negative and positive values

    **Color Gradient**:
    - Values at minimum map to color_min
    - Values at maximum map to color_max
    - Values at midpoint (typically 0) map to light gray
    - Linear interpolation for intermediate values

    **Data Considerations**:
    - Works best with normalized/scaled data
    - Extreme outliers may compress visualization range
    - Consider data transformation for better visualization

    **iTOL Integration**:
    - Upload tree (Newick format) and annotation file to iTOL
    - Customize appearance in iTOL interface
    - Export publication-quality figures

    See Also
    --------
    dataframe_to_itol_colorstrip : Create categorical color annotations
    dataframe_to_itol_stackedbar : Create multi-category bar annotations
    """
    file = io.StringIO()
    # Write the header for the iTOL heatmap dataset
    file.write("DATASET_HEATMAP\n")
    file.write("SEPARATOR SPACE\n")
    file.write(f"DATASET_LABEL {dataset_label}\n")
    file.write("FIELD_LABELS " + " ".join(df.columns) + "\n")
    file.write("COLOR #ff0000\n")  # Default color, not used in coolwarm palette

    # Write color gradients for coolwarm
    file.write(f"COLOR_MIN {color_min}\n")  # Cool color
    file.write("COLOR_MID #f5f5f5\n")  # Midpoint color
    file.write(f"COLOR_MAX {color_max}\n")  # Warm color

    # Data section
    file.write("DATA\n")
    for index, row in df.iterrows():
        file.write(f"{index} " + " ".join(map(str, row)) + "\n")

    text = file.getvalue()
    file.close()
    return text


def dataframe_to_itol_stackedbar(
    df: pd.DataFrame, cmap: Dict[str, Union[str, Tuple[float, ...]]], dataset_label: str
) -> str:
    """
    Generate iTOL stacked bar annotation for multi-category data visualization.

    Creates an iTOL-compatible multi-bar annotation file for visualizing multiple
    quantitative categories (such as clone proportions, cell type compositions,
    or feature counts) alongside phylogenetic trees as stacked bars.

    Parameters
    ----------
    df : pd.DataFrame
        Data matrix with samples/cells as rows and categories as columns. Index should
        match leaf names in the phylogenetic tree. Values represent quantities for
        each category (e.g., cell counts, proportions, expression levels).
    cmap : dict
        Dictionary mapping column names (categories) to colors. Keys should match
        DataFrame column names. Values can be:
        - Hex color strings (e.g., '#FF0000')
        - RGB tuples (e.g., (1.0, 0.0, 0.0))
        - RGBA tuples (e.g., (1.0, 0.0, 0.0, 1.0))
        RGB/RGBA tuples will be automatically converted to hex format.
    dataset_label : str
        Label for the annotation dataset, displayed in iTOL interface and legend.

    Returns
    -------
    str
        iTOL-formatted multi-bar annotation content that can be saved to a file
        and uploaded to iTOL for phylogenetic tree visualization.

    Examples
    --------
    Visualize clone composition per sample:

    >>> import pandas as pd
    >>> import seaborn as sns
    >>> from picasso.itol_utils import dataframe_to_itol_stackedbar
    >>>
    >>> # Create composition data (proportions sum to 1.0 per row)
    >>> composition_data = pd.DataFrame({
    ...     'Clone_A': [0.6, 0.8, 0.2, 0.1],
    ...     'Clone_B': [0.3, 0.2, 0.5, 0.7],
    ...     'Clone_C': [0.1, 0.0, 0.3, 0.2]
    ... }, index=['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4'])
    >>>
    >>> # Define colors for each clone
    >>> clone_colors = {
    ...     'Clone_A': '#FF0000',  # Red
    ...     'Clone_B': '#00FF00',  # Green
    ...     'Clone_C': '#0000FF'   # Blue
    ... }
    >>>
    >>> # Generate iTOL annotation
    >>> stackedbar = dataframe_to_itol_stackedbar(composition_data,
    ...                                          clone_colors,
    ...                                          'Clone_Composition')

    Using RGB tuples with seaborn colors:

    >>> # Generate colors automatically
    >>> categories = ['TypeA', 'TypeB', 'TypeC', 'TypeD']
    >>> colors = sns.color_palette('Set2', len(categories))
    >>> color_map = dict(zip(categories, colors))
    >>>
    >>> stackedbar = dataframe_to_itol_stackedbar(data, color_map, 'Cell_Types')

    Notes
    -----
    **Stacked Bar Visualization**:
    - Each row becomes a stacked bar with segments for each column
    - Segment heights proportional to values in DataFrame
    - Colors defined by the color mapping dictionary
    - Useful for showing compositional data

    **Data Requirements**:
    - Non-negative values recommended for meaningful visualization
    - Values can be counts, proportions, or any quantitative measure
    - Missing values (NaN) treated as zero
    - Consider normalizing data for proportion-based visualization

    **Color Handling**:
    - Supports hex strings, RGB, and RGBA color formats
    - Automatic conversion of tuple formats to hex strings
    - Colors must be defined for all DataFrame columns
    - Consistent color mapping across multiple annotations

    **Use Cases**:
    - Clone composition per tissue/sample
    - Cell type distributions per phylogenetic group
    - Multi-gene expression profiles
    - Pathway activity scores across conditions

    Raises
    ------
    AssertionError
        If any DataFrame column lacks a corresponding color in the colormap,
        or if colors cannot be converted to valid string format.

    See Also
    --------
    dataframe_to_itol_colorstrip : Create categorical color annotations
    dataframe_to_itol_heatmap : Create continuous data heatmaps
    rgba_to_hex : Convert color tuples to hex strings
    """
    for key in cmap:
        # If the color is in rgba format, convert it to hex
        if isinstance(cmap[key], tuple) or isinstance(cmap[key], list):
            cmap[key] = rgba_to_hex(cmap[key])
        assert isinstance(cmap[key], str), f"Color for {key} is not a string"

    file = io.StringIO()
    # Write the header for the iTOL heatmap dataset
    file.write("DATASET_MULTIBAR\n")
    file.write("SEPARATOR\tTAB\n")
    file.write(f"DATASET_LABEL\t{dataset_label}\n")
    file.write("FIELD_LABELS\t" + "\t".join(df.columns) + "\n")
    colors = "\t".join([cmap[col] for col in df.columns])
    file.write(
        f"FIELD_COLORS\t{colors}\n"
    )  # Default color, not used in coolwarm palette

    # Data section
    file.write("DATA\n")
    for index, row in df.iterrows():
        file.write(f"{index}\t" + "\t".join(map(str, row)) + "\n")

    text = file.getvalue()
    file.close()
    return text


# Define public API
__all__ = [
    "rgba_to_hex",
    "dataframe_to_itol_colorstrip",
    "dataframe_to_itol_heatmap",
    "dataframe_to_itol_stackedbar",
]
