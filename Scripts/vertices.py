import os
import re
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy.spatial.distance import cdist

def extract_timestamp(filename):
    # Extract timestamp from filename (assuming format: 2d_grain30000.dat)
    print(filename)
    match = re.search(r'2d_grain(\d+)\.dat', filename)
    if match:
        return int(match.group(1))
    return 0

def identify_vertices(grain_grid, max_x, max_y, vulnerable_zones):
    vertices = []
    
    # Identify grain boundaries and vertices
    for x in range(max_x):
        for y in range(max_y):
            if grain_grid[x, y] == -1:
                continue
            
            # Check for grain boundaries
            is_boundary = False
            grain_id = grain_grid[x, y]
            
            neighbors = [
                (x + dx, y + dy)
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]
                if 0 <= x + dx < max_x and 0 <= y + dy < max_y
            ]
            
            for nx, ny in neighbors:
                if grain_grid[nx, ny] != grain_id and grain_grid[nx, ny] != -1:
                    is_boundary = True
                    break
            
            # If this is a boundary pixel, check if it's a vertex
            if is_boundary:
                adjacent_grains = set()
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        nx, ny = x + dx, y + dy
                        if 0 <= nx < max_x and 0 <= ny < max_y and grain_grid[nx, ny] != -1:
                            adjacent_grains.add(grain_grid[nx, ny])
                
                # If the pixel touches 3 or more grains, it's a vertex
                if len(adjacent_grains) >= 3:
                    is_vulnerable = vulnerable_zones[x, y] if vulnerable_zones is not None else False
                    vertices.append((x, y, is_vulnerable, tuple(sorted(adjacent_grains))))
    
    return vertices

def match_vertices_between_frames(prev_vertices, curr_vertices, max_distance=10):
    """Match vertices between consecutive frames based on proximity and grain connections"""
    if not prev_vertices:
        return {}, {}
    
    # Extract coordinates
    prev_coords = np.array([(x, y) for x, y, _, _ in prev_vertices])
    curr_coords = np.array([(x, y) for x, y, _, _ in curr_vertices])
    
    # Create dictionaries for faster lookup
    prev_dict = {(x, y): (is_vuln, grains) for x, y, is_vuln, grains in prev_vertices}
    curr_dict = {(x, y): (is_vuln, grains) for x, y, is_vuln, grains in curr_vertices}
    
    # Calculate distances between all pairs of vertices
    distances = cdist(prev_coords, curr_coords)
    
    # Track matched vertices
    matches = {}
    reverse_matches = {}
    
    # For each previous vertex, find the closest current vertex within threshold
    for i, (px, py, _, prev_grains) in enumerate(prev_vertices):
        sorted_indices = np.argsort(distances[i])
        
        for j in sorted_indices:
            if distances[i][j] <= max_distance:
                cx, cy = curr_coords[j]
                _, curr_grains = curr_dict[(cx, cy)]
                
                # Check if there's some grain overlap to confirm it's the same vertex
                if len(set(prev_grains).intersection(set(curr_grains))) >= 1:  # Relaxed condition (at least 1 common grain)
                    matches[(px, py)] = (cx, cy)
                    reverse_matches[(cx, cy)] = (px, py)
                    break
    
    return matches, reverse_matches

def track_vertices_and_calculate_velocities(input_directory):
    print("Tracking vertices and calculating velocities...")
    
    # Collect all dat files with timestamps
    all_files = []
    for folder, _, files in os.walk(input_directory):
        for file in files:
            if file.startswith("2d_grain") and file.endswith(".dat"):
                timestamp = extract_timestamp(file)
                if timestamp > 0:  # Valid timestamp found
                    all_files.append((os.path.join(folder, file), timestamp))
    
    if not all_files:
        print(f"No valid timestamped files found in {input_directory}")
        print("Example filename expected: 2d_grain30000.dat")
        return {}
    
    # Sort files by timestamp
    all_files.sort(key=lambda x: x[1])
    print(f"Found {len(all_files)} timestamped files. Processing...")
    
    # Track vertices over time
    vertex_positions = defaultdict(list)  # vertex_id -> [(timestamp, x, y, is_vulnerable), ...]
    next_vertex_id = 0
    
    prev_vertices = []
    prev_timestamp = None
    vertex_to_id_map = {}  # Map (x, y) coordinates to vertex ID
    
    for file_path, timestamp in all_files:
        print(f"Processing: {os.path.basename(file_path)} (timestamp: {timestamp})")
        
        try:
            # Load data
            data = []
            with open(file_path, "r") as infile:
                for line in infile:
                    columns = line.strip().split()
                    try:
                        x = int(columns[0])
                        y = int(columns[1])
                        grain_id = int(float(columns[-4]))
                        value = int(columns[-1])
                        data.append((x, y, grain_id, value))
                    except (ValueError, IndexError):
                        continue
            
            if not data:
                print(f"No valid data found in {file_path}")
                continue
                
            try:
                max_x = max(row[0] for row in data) + 1
                max_y = max(row[1] for row in data) + 1
            except ValueError:
                print(f"Error calculating grid dimensions for {file_path}")
                continue
            
            # Create grids
            grain_grid = np.full((max_x, max_y), -1, dtype=int)
            vulnerable_zones = np.zeros((max_x, max_y), dtype=bool)
            
            # Fill grids
            for x, y, grain_id, value in data:
                grain_grid[x, y] = grain_id
                
                if value == 1:
                    for i in range(-2, 3):
                        for j in range(-2, 3):
                            if 0 <= x + i < max_x and 0 <= y + j < max_y:
                                vulnerable_zones[x + i, y + j] = True
            
            # Identify vertices in current frame
            curr_vertices = identify_vertices(grain_grid, max_x, max_y, vulnerable_zones)
            
            if not curr_vertices:
                print(f"No vertices found in {file_path}")
                continue
                
            print(f"  Found {len(curr_vertices)} vertices")
            
            # Match vertices with previous frame
            vertex_matches, reverse_matches = match_vertices_between_frames(prev_vertices, curr_vertices)
            
            # Update vertex tracking
            for cx, cy, is_vuln, grains in curr_vertices:
                if (cx, cy) in reverse_matches:
                    # This vertex was matched with a previous one
                    px, py = reverse_matches[(cx, cy)]
                    
                    if (px, py) in vertex_to_id_map:
                        vertex_id = vertex_to_id_map[(px, py)]
                        vertex_positions[vertex_id].append((timestamp, cx, cy, is_vuln))
                        vertex_to_id_map[(cx, cy)] = vertex_id  # Update map with new position
                else:
                    # This is a new vertex
                    vertex_id = next_vertex_id
                    next_vertex_id += 1
                    vertex_positions[vertex_id].append((timestamp, cx, cy, is_vuln))
                    vertex_to_id_map[(cx, cy)] = vertex_id
            
            # Save current vertices for next iteration
            prev_vertices = curr_vertices
            prev_timestamp = timestamp
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    print(f"Tracking complete. Found {len(vertex_positions)} vertices across all timestamps.")
    
    # Calculate velocities
    vertex_velocities = {}
    
    for vertex_id, positions in vertex_positions.items():
        if len(positions) >= 2:  # Need at least 2 positions to calculate velocity
            velocities = []
            for i in range(1, len(positions)):
                prev_time, prev_x, prev_y, _ = positions[i-1]
                curr_time, curr_x, curr_y, is_vuln = positions[i]
                
                # Skip if timestamps are the same to avoid division by zero
                if curr_time == prev_time:
                    continue
                
                # Calculate velocity (pixels per time unit)
                distance = np.sqrt((curr_x - prev_x)**2 + (curr_y - prev_y)**2)
                time_diff = curr_time - prev_time
                velocity = distance / time_diff
                
                velocities.append((curr_time, velocity, is_vuln))
            
            if velocities:
                vertex_velocities[vertex_id] = velocities
    
    print(f"Velocity calculation complete. {len(vertex_velocities)} vertices have valid velocity data.")
    return vertex_velocities

def plot_vertex_velocities(vertex_velocities, output_directory):
    """Create plots showing velocity profiles of all vertices"""
    if not vertex_velocities:
        print("No vertex velocity data to plot")
        return
    
    print(f"Creating velocity plots for {len(vertex_velocities)} vertices...")
    
    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)
    
    # Create a plot for all vertices
    plt.figure(figsize=(12, 8))
    
    # Split into vulnerable and non-vulnerable vertices for better visualization
    vulnerable_data = []
    normal_data = []
    
    for vertex_id, velocities in vertex_velocities.items():
        times = [v[0] for v in velocities]
        speeds = [v[1] for v in velocities]
        is_vulnerable = any(v[2] for v in velocities)
        
        if is_vulnerable:
            vulnerable_data.append((times, speeds, vertex_id))
        else:
            normal_data.append((times, speeds, vertex_id))
    
    # Plot normal vertices
    for times, speeds, vertex_id in normal_data:
        plt.plot(times, speeds, 'b-', alpha=0.3, linewidth=0.8)
    
    # Plot vulnerable vertices
    for times, speeds, vertex_id in vulnerable_data:
        plt.plot(times, speeds, 'r-', alpha=0.5, linewidth=1.2)
    
    # Add legend
    plt.plot([], [], 'b-', label='Normal Vertices')
    plt.plot([], [], 'r-', label='Vulnerable Zone Vertices')
    
    plt.xlabel('Timestamp')
    plt.ylabel('Velocity (pixels/time unit)')
    plt.title('Vertex Velocity Profiles Over Time')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Save the plot
    velocity_profile_path = os.path.join(output_directory, 'vertex_velocity_profiles.png')
    plt.savefig(velocity_profile_path, dpi=300, bbox_inches='tight')
    print(f"Saved velocity profile plot to {velocity_profile_path}")
    
    # Collect all velocities for statistics and comparison
    vulnerable_velocities = []
    normal_velocities = []
    
    for vertex_id, velocities in vertex_velocities.items():
        for _, speed, is_vuln in velocities:
            if is_vuln:
                vulnerable_velocities.append(speed)
            else:
                normal_velocities.append(speed)
    
    # Create boxplot comparing normal vs vulnerable vertices
    plt.figure(figsize=(10, 8))
    
    # Create boxplot
    box_data = [normal_velocities, vulnerable_velocities]
    plt.boxplot(box_data, labels=['Normal Vertices', 'Vulnerable Zone Vertices'])
    plt.ylabel('Velocity Distribution (pixels/time unit)')
    plt.title('Comparison of Vertex Velocities')
    plt.grid(True, alpha=0.3)
    
    # Save the comparison plot
    comparison_plot_path = os.path.join(output_directory, 'vertex_velocity_comparison.png')
    plt.savefig(comparison_plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved velocity comparison plot to {comparison_plot_path}")
    
    # Histogram of velocities
    plt.figure(figsize=(12, 8))
    
    max_velocity = max(max(vulnerable_velocities, default=0), max(normal_velocities, default=0))
    bins = np.linspace(0, max_velocity * 1.1 if max_velocity > 0 else 10, 30)
    
    plt.hist(normal_velocities, bins=bins, alpha=0.5, label='Normal Vertices', color='blue')
    plt.hist(vulnerable_velocities, bins=bins, alpha=0.5, label='Vulnerable Zone Vertices', color='red')
    
    plt.xlabel('Velocity (pixels/time unit)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Vertex Velocities')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save the histogram
    histogram_path = os.path.join(output_directory, 'vertex_velocity_histogram.png')
    plt.savefig(histogram_path, dpi=300, bbox_inches='tight')
    print(f"Saved velocity histogram to {histogram_path}")
    
    # Create summary statistics
    print("\nVelocity Statistics:")
    print(f"Normal vertices count: {len(normal_velocities)}")
    print(f"Vulnerable vertices count: {len(vulnerable_velocities)}")
    
    if normal_velocities:
        print(f"Normal vertices - Mean velocity: {np.mean(normal_velocities):.4f}, Median: {np.median(normal_velocities):.4f}, Max: {np.max(normal_velocities):.4f}")
    
    if vulnerable_velocities:
        print(f"Vulnerable vertices - Mean velocity: {np.mean(vulnerable_velocities):.4f}, Median: {np.median(vulnerable_velocities):.4f}, Max: {np.max(vulnerable_velocities):.4f}")

    # Save statistics to a file
    stats_path = os.path.join(output_directory, 'velocity_statistics.txt')
    with open(stats_path, 'w') as stats_file:
        stats_file.write("Vertex Velocity Statistics\n")
        stats_file.write("==========================\n\n")
        stats_file.write(f"Total vertices tracked: {len(vertex_velocities)}\n")
        stats_file.write(f"Normal vertices count: {len(normal_velocities)}\n")
        stats_file.write(f"Vulnerable vertices count: {len(vulnerable_velocities)}\n\n")
        
        if normal_velocities:
            stats_file.write("Normal vertices statistics:\n")
            stats_file.write(f"  Mean velocity: {np.mean(normal_velocities):.4f}\n")
            stats_file.write(f"  Median velocity: {np.median(normal_velocities):.4f}\n")
            stats_file.write(f"  Max velocity: {np.max(normal_velocities):.4f}\n")
            stats_file.write(f"  Min velocity: {np.min(normal_velocities):.4f}\n")
            stats_file.write(f"  Standard deviation: {np.std(normal_velocities):.4f}\n\n")
        
        if vulnerable_velocities:
            stats_file.write("Vulnerable vertices statistics:\n")
            stats_file.write(f"  Mean velocity: {np.mean(vulnerable_velocities):.4f}\n")
            stats_file.write(f"  Median velocity: {np.median(vulnerable_velocities):.4f}\n")
            stats_file.write(f"  Max velocity: {np.max(vulnerable_velocities):.4f}\n")
            stats_file.write(f"  Min velocity: {np.min(vulnerable_velocities):.4f}\n")
            stats_file.write(f"  Standard deviation: {np.std(vulnerable_velocities):.4f}\n")
    
    print(f"Saved detailed statistics to {stats_path}")

def main():
    # Update paths
    input_directory = r"C:\project\Btp\Dewetting-of-thin-films\processed\Dewetting-of-thin-films\data\ms_100_24\output\grain"
    output_directory = "./grain_vertices_velocity"
    
    # Track vertices and calculate velocities
    vertex_velocities = track_vertices_and_calculate_velocities(input_directory)
    
    # Plot velocity profiles
    plot_vertex_velocities(vertex_velocities, output_directory)

if __name__ == "__main__":
    main()