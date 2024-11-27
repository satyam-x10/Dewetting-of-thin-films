import pandas as pd

def calculate_interface_depth(file_path, output_path, surface_z=39, mod_threshold=0.5):
    """
    Calculate depth based on the vapor-film interface.

    Parameters:
        file_path (str): Path to the input data file.
        output_path (str): Path to save the updated data with depth values.
        surface_z (float): Reference surface level (z-coordinate).
        mod_threshold (float): Threshold for identifying the vapor-film interface.

    Returns:
        pd.DataFrame: Updated DataFrame with interface depth for each (x, y).
    """
    # Load the data
    print(f"Loading data from: {file_path}")
    data = pd.read_csv(file_path, delim_whitespace=True, header=None, 
                       names=['x', 'y', 'z', 'mod', 'grain_id', 'sum_eta2', 'sum_eta', 'vulnerable'])
    print(f"Data loaded. Total rows: {len(data)}")

    # Verify z column range
    z_min, z_max = data['z'].min(), data['z'].max()
    print(f"z column range: min={z_min}, max={z_max}")

    # Group data by (x, y) to analyze transitions in z for each (x, y)
    depths = []
    for (x, y), group in data.groupby(['x', 'y']):
        print(f"\nProcessing (x={x}, y={y}) with {len(group)} points.")
        
        # Sort by z to ensure correct order
        group = group.sort_values(by='z')
        print(f"Sorted group for (x={x}, y={y}):")
        print(group[['z', 'mod']].head())  # Display the first few rows for debugging

        # Find the interface where mod transitions from >0.5 to <=0.5
        interface_row = None
        for i in range(1, len(group)):
            if group.iloc[i-1]['mod'] > mod_threshold and group.iloc[i]['mod'] <= mod_threshold:
                interface_row = group.iloc[i]
                print(f"  Interface found at z={interface_row['z']} for (x={x}, y={y}).")
                break
        
        # Calculate depth if an interface is found
        if interface_row is not None:
            z_interface = interface_row['z']
            depth = z_interface - surface_z
            depths.append({'x': x, 'y': y, 'z_interface': z_interface, 'depth': depth})
        else:
            print(f"  No interface found for (x={x}, y={y}).")
            depths.append({'x': x, 'y': y, 'z_interface': None, 'depth': None})

    # Create a DataFrame for the depths
    depth_df = pd.DataFrame(depths)

    # Save the depths to a file
    depth_df.to_csv(output_path, index=False, sep='\t')
    print(f"Depth data saved to: {output_path}")

    return depth_df

# Input file path (uploaded file path)
input_file = r'G:\BTP\Dewetting-of-thin-films\processed\24h_1000\output\grain\3d_grain108000.dat'
output_file = r'G:\BTP\Dewetting-of-thin-films\depths\3d_grain108000_depths.dat'  # Path to save the depth data

# Calculate depths based on the vapor-film interface
depth_data = calculate_interface_depth(input_file, output_file, surface_z=39, mod_threshold=0.5)
