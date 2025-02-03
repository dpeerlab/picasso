import io

def rgba_to_hex(rgba):
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
    red_hex = format(int(red * 255), '02X')
    green_hex = format(int(green * 255), '02X')
    blue_hex = format(int(blue * 255), '02X')

    # Concatenate the hexadecimal values
    hex_color = f"#{red_hex}{green_hex}{blue_hex}"

    return hex_color

def dataframe_to_itol_colorstrip(series, cmap, dataset_label):
    """
    Convert a pandas DataFrame into an iTOL colorstrip annotation file.
    :param series: Pandas Series with index as leaf labels and values as data points.
    :param cmap: Dictionary mapping data points to colors.
    """

    for key in cmap:
        # If the color is in rgba format, convert it to hex
        if isinstance(cmap[key], tuple) or isinstance(cmap[key], list):
            cmap[key] = rgba_to_hex(cmap[key])
        assert isinstance(cmap[key], str), f"Color for {key} is not a string"

    # Create the annotations file and write it to buffer

    f = io.StringIO()
    f.write('DATASET_COLORSTRIP\n')
    f.write('SEPARATOR TAB\n')
    f.write(f'DATASET_LABEL\t{dataset_label}\n')
    f.write('COLOR\t#ff0000\n')
    f.write(f'LEGEND_TITLE\t{dataset_label}\n')
    f.write('LEGEND_SHAPES\t1\n')
    f.write('LEGEND_COLORS\t#ff0000\n')
    f.write(f'LEGEND_LABELS\t{dataset_label}\n')
    f.write('DATA\n')
    for leaf in series.index:
        lineage = series.loc[leaf]
        f.write(f'{leaf}\t{cmap[lineage]}\t{lineage}\n')
    text = f.getvalue()
    f.close()
    return text


def dataframe_to_itol_heatmap(df, dataset_label="CNVs", color_min='#3f4c8a', color_max='#b40426'):
    """
    Convert a pandas DataFrame into an iTOL heatmap annotation file.
    :param df: Pandas DataFrame with rows as leaf labels and columns as data points.
    :param dataset_label: Label for the dataset.
    :param color_min: Color for the minimum value (as a hex string).
    :param color_max: Color for the maximum value (as a hex string).
    :param output_file: Path to the output file.
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

def dataframe_to_itol_stackedbar(df, cmap, dataset_label):
    """
    Convert a pandas DataFrame into an iTOL heatmap annotation file.

    :param df: Pandas DataFrame with rows as leaf labels and columns as data points.
    :param output_file: Path to the output file.
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
    file.write(f"FIELD_COLORS\t{colors}\n")  # Default color, not used in coolwarm palette

    # Data section
    file.write("DATA\n")
    for index, row in df.iterrows():
        file.write(f"{index}\t" + "\t".join(map(str, row)) + "\n")

    text = file.getvalue()
    file.close()
    return text



